#------------------------------------------------------------------------------------------------------------
#' Transforms an assay in an \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @inheritParams generic_scMethrix_function
#' @param h5_temp string; temporary directory to store hdf5
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' transform_assay(scMethrix_data,assay="score",new_assay="plus1",trans=function(x){x+1})
#' @export
transform_assay <- function(scm, assay = NULL, new_assay = NULL, trans = NULL, h5_temp = NULL) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (typeof(trans) != "closure") {
    stop("A valid transform function must be specified.", call. = FALSE)
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }
  
  if (new_assay %in% SummarizedExperiment::assayNames(scm)) {
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)
  }
  
  if (is_h5(scm)) {
    
    if (is.null(h5_temp)) {h5_temp <- tempdir()}
    
    grid <- DelayedArray::RegularArrayGrid(refdim = dim(scm),
                                           spacings = c(length(scm), 1L)) 
    
    trans_sink <- HDF5Array::HDF5RealizationSink(dim = dim(scm),
                                                 dimnames = list(NULL,row.names(colData(scm))), type = "integer",
                                                 filepath = tempfile(pattern="trans_sink_",tmpdir=h5_temp),
                                                 name = new_assay, level = 6)
    
    blocs <- DelayedArray::blockApply(get_matrix(scm,assay=assay), grid = grid, FUN = trans)
    
    for(i in 1:length(blocs)) {
      DelayedArray::write_block(block = as.matrix(blocs[[i]]), viewport = grid[[i]], sink = trans_sink)
    }
    
    rm(blocs)
    s <- as(trans_sink, "HDF5Matrix")
    
  } else {
    s <- as.data.table(get_matrix(scm,assay=assay))
    s[, names(s) := lapply(.SD, trans)]
    s <- as(s, "matrix")
  }
  
  assays(scm)[[new_assay]] <- s
  
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Bins the ranges of an \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @inheritParams generic_scMethrix_function
#' @param bin_size The size of each bin. First bin will begin at the start position of the first genomic
#' region on the chromosome
#' @param trans The transforms for each assay. Must be a named vector of functions (closure). 
#' Default = mean(.., na.rm=TRUE)
#' @param h5_dir directory to store H5 based object
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' trans <- c(score = function(x) mean(x,na.rm=TRUE),counts = function(x) sum(x,na.rm=TRUE))
#' bin_scMethrix(scMethrix_data,trans = trans)
#' @export
bin_scMethrix <- function(scm, bin_size = 100000, trans = NULL, h5_dir = NULL) {
  
  if (is.null(trans)) {
    trans <- c(score = function(x) mean(x,na.rm=TRUE),counts = function(x) sum(x,na.rm=TRUE))
  }
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is_h5(scm) && is.null(h5_dir)) stop("Output directory must be specified", call. = FALSE)
  
  bins <- bin_granges(rowRanges(scm),bin_size = bin_size)
  bins <- subsetByOverlaps(bins,rowRanges(scm))
  
  sites <- sapply(1:length(bins),function (i) {
    (findOverlaps(rowRanges(scm),bins[i]))@from
  })
  
  assays <- list()
  
  for (n in 1:length(assays(scm))) {
    
    name <- SummarizedExperiment::assayNames(scm)[n]
    
    tryCatch(
      expr = {op <- trans[[name]]},
      error = function(e){op <- function(x) mean(x,na.rm=TRUE)}
    )
    
    vals <- lapply(sites,function (i) {
      apply(get_matrix(scm[i,],assay=name),2,op)
    })
    
    setDT(vals, key=names(vals[[1]]))
    vals <- data.table::transpose(vals)
    colnames(vals) <- rownames(colData(scm))
    
    assays[[name]] <- as(vals,class(get_matrix(scm,assay=name)))
  }
  
  if (is_h5(scm)) {
    
    m_obj <- create_scMethrix(assays = assays, rowRanges=bins, is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),
                              replace = replace)
    
  } else {
    
    m_obj <- create_scMethrix(assays = assays, rowRanges=bins, is_hdf5 = FALSE, 
                              genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),)
    
  }
}

#------------------------------------------------------------------------------------------------------------
#' Imputes the NA values of a \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @param threshold The value for cutoff in the "score" assay to determine methylated or unmethylated status. 
#' Default = 50
#' @inheritParams generic_scMethrix_function
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' @export
#' @import Melissa
impute_by_melissa <- function (scm, threshold = 50, assay = "score", new_assay = "impute") {
  
  . <- NULL
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  scm <- transform_assay(scm, assay = assay, new_assay = "binary", trans = binarize)

  # Convert Granges to genomic interval [-1,1]
  chrom_size <- scm@metadata$chrom_size
  chrom_start <- sapply(coverage(scm), function(x) {x@lengths[1]}+1)
  interval <- as.data.table(rowRanges(scm))[,.(seqnames,start)]
  
  for (i in 1: length(chrom_size)) {
    # The final square bracket at line below is added due to bug:
    # https://stackoverflow.com/questions/34667536/have-to-call-variable-twice-before-evaluated
    interval[seqnames==names(chrom_start)[i], interval := round((2*((start-chrom_start[i])/chrom_size[i]))-1, digits = 7)][]
  }
  
  mcols(scm) <- cbind(mcols(scm),interval = interval[,interval])
  
  # Put into Melissa met format
  loc <- as.matrix(rowRanges(scm)@seqnames)
  loc <- matrix(loc,dimnames = list(paste0(rowRanges(scm)@seqnames,":",interval$start)))
  cells <- list()
  
  for (n in 1:ncol(scm)) {
    
    met <- matrix(c(mcols(scm)$interval,get_matrix(scm[,n],assay="binary")),
                  ncol = 2,dimnames = list(rownames(loc),NULL))
    mets <- lapply(rowRanges(scm)@seqnames@values, function(x) met[which(loc == x),,drop=FALSE])
    mets <- lapply(mets, function (m) m[-which(m[,2]==-1),,drop=FALSE])
    
    cells[[colnames(scm)[n]]] <- mets
  }
  
  # # Parameter options
  opts <- list()
  # # Load cell filenames
  # opts$met_files <- NULL
  # opts$cell_id <- colnames(scm)
  # opts$is_centre  <- is_centre   # Whether genomic region is already pre-centred
  # opts$is_window  <- is_window   # Use predefined window region
  # opts$upstream   <- upstream    # Upstream of centre
  # opts$downstream <- downstream  # Downstream of centre
  # opts$chrom_size <- chrom_size_file  # Chromosome size file
  # opts$chr_discarded <- chr_discarded # Chromosomes to discard
  # opts$cov        <- cov         # Regions with at least n CpGs
  # opts$sd_thresh  <- sd_thresh   # Variance of methylation within region
  # 
  melissa_obj <- structure(list(met = cells, anno_region = NULL, opts = opts),
                           class = "melissa_data_obj")
  
  # Do the imputation
  basis_obj <- BPRMeth::create_rbf_object(M = 4)
  test_obj <- Melissa::partition_dataset(melissa_obj)
  
  set.seed(123)
  melissa_obj <- Melissa::melissa(X = melissa_obj$met, basis = basis_obj,K = min(ncol(scm)-1,2))#,
                                  #vb_max_iter = 30, vb_init_nstart = 1, 
                                  #is_parallel = FALSE)
  # plot_melissa_profiles(melissa_obj = melissa_obj, region = 1, 
  #                       title = "Methylation profiles for region 25")
  
  
  
  imputation_obj <- Melissa::impute_test_met(obj = melissa_obj, test = test_obj)
                                    
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Imputes the an \code{\link{scMethrix}}object using a iterative PCA model
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object.
#' @param n_pc integer > 0; The number of principal components to use. This can be a range - e.g. c(1,5) - or a single value. 
#' For a range, the optimal value will be estimated; this is time-intensive.
#' @param ... Additional arguments for missMDA::imputePCA
#' @inheritParams transform_assay
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' impute_by_iPCA(scMethrix_data, assay = "score", new_assay = "impute", n_pc = 2)
#' @export
#' @references Bro, R., Kjeldahl, K. Smilde, A. K. and Kiers, H. A. L. (2008) Cross-validation of component models: A critical look at current methods. Analytical and Bioanalytical Chemistry, 5, 1241-1251.
#' @references Josse, J. and Husson, F. (2011). Selecting the number of components in PCA using cross-validation approximations. Computational Statistics and Data Analysis. 56 (6), pp. 1869-1879.
impute_by_iPCA <- function(scm = NULL, assay = "score", new_assay = "impute", n_pc = 2, ...) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }
  
  if (new_assay %in% SummarizedExperiment::assayNames(scm)) {
    if (new_assay == "score") stop("Cannot overwrite the score assay")
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)
  }
  
  if (is_h5(scm)) {
    stop("HDF5-based objects cannot be imputed by iPCA", call. = FALSE)
  }
  
  if (length(n_pc) > 1) {
    warning("Caution: n_pc is given as range. This can be very time-intensive.")
    n_pc <- missMDA::estim_ncpPCA(get_matrix(scm,assay = assay),ncp.min = n_pc[1], ncp.max = n_pc[2], 
                         method.cv = "Kfold", verbose = TRUE)
    n_pc <- n_pc$ncp
  }
  
  impute <- missMDA::imputePCA(get_matrix(scm,assay = assay), ncp = n_pc, ...)
  assays(scm)[[new_assay]] <- as(impute$completeObs,class(get_matrix(scm,assay=assay))[[1]])
  
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Imputes the an \code{\link{scMethrix}}object using a random forest model
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object.
#' @param ... Additional arguments for missForest::missForest
#' @inheritParams transform_assay
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' impute_by_RF(scMethrix_data, assay = "score", new_assay = "impute")
#' @export
#' @import missForest
impute_by_RF <- function(scm = NULL, assay = "score", new_assay = "impute", ...) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }
  
  if (new_assay %in% SummarizedExperiment::assayNames(scm)) {
    if (new_assay == "score") stop("Cannot overwrite the score assay")
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)
  }
  
  if (is_h5(scm)) {
    stop("HDF5-based objects cannot be imputed by RF", call. = FALSE)
  }
  
  impute <- missForest::missForest(get_matrix(scm,assay = assay))#, ...)
  assays(scm)[[new_assay]] <- as(impute$ximp,class(get_matrix(scm,assay=assay))[[1]])
  
  return(scm)
}


#------------------------------------------------------------------------------------------------------------
#' Imputes the an \code{\link{scMethrix}}object using a random forest model
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object.
#' @param ... Additional arguments for missForest::missForest
#' @inheritParams impute::impute.knn
#' @inheritParams transform_assay
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' impute_by_RF(scMethrix_data, assay = "score", new_assay = "impute")
#' @export
#' @import impute
impute_by_kNN <- function(scm = NULL, assay = "score", new_assay = "impute", k = 10, ...) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }
  
  if (new_assay %in% SummarizedExperiment::assayNames(scm)) {
    if (new_assay == "score") stop("Cannot overwrite the score assay")
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)
  }
  
  if (is_h5(scm)) {
    stop("HDF5-based objects cannot be imputed by RF", call. = FALSE)
  }
  
  impute <- impute::impute.knn(get_matrix(scm,assay = assay), k = min(k,ncol(scm)), 
                       rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=123)
  assays(scm)[[new_assay]] <- as(impute$data,class(get_matrix(scm,assay=assay))[[1]])
  
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Splits an experiment into two for use as a training and test set
#' @details Does stuff
#' @param training_prop numeric; The size of the training set as a proportion of the experiment (0 to 1)
#' For a range, the optimal value will be estimated; this is time-intensive.
#' @inheritParams generic_scMethrix_function
#' @return list; two \code{\link{scMethrix}} objects names 'training' and 'test'
#' @examples
#' data('scMethrix_data')
#' generate_training_set(scMethrix_data, training_prop = 0.2)
#' @export
generate_training_set <- function(scm = NULL, training_prop = 0.2) {
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (training_prop > 1 || training_prop < 0) {
    stop("training_prop must in the range of [0,1]", call. = FALSE)
  }
  set.seed(123)
  idx <- sort(sample(1:nrow(scm),floor(nrow(scm)*training_prop)))
  
  training <- scm[idx,]
  test <- scm[setdiff(1:nrow(scm),idx),]
  
  return(list(training = training,test = test))
}