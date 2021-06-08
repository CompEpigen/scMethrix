#------------------------------------------------------------------------------------------------------------
#' Transforms an assay in an \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @param scm A \code{\link{scMethrix}} object
#' @param assay String name of an existing assay
#' @param new_assay String name of transformed assay
#' @param trans The transformation function
#' @param h5_temp temporary directory to store hdf5
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
#' @param scm A \code{\link{scMethrix}} object
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
  
  assys <- list()
  
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
    
    assys[[name]] <- as(vals,class(get_matrix(scm,assay=name)))
  }
  
  if (is_h5(scm)) {
    
    m_obj <- create_scMethrix(assays = assys, rowRanges=bins, is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),
                              replace = replace)
    
  } else {
    
    m_obj <- create_scMethrix(assays = assys, rowRanges=bins, is_hdf5 = FALSE, 
                              genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),)
    
  }
}

#------------------------------------------------------------------------------------------------------------
#' Imputes the NA values of a \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @param scm A \code{\link{scMethrix}} object
#' @param threshold The value for cutoff in the "score" assay to determine methylated or unmethylated status. 
#' Default = 50
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' @export
impute_by_melissa <- function (scm, threshold = 50) {
  
  . <- NULL
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  scm <- transform_assay(scm, assay = "score", new_assay = "binary",trans = binarize)
  
  melissa_obj <- generate_melissa_object(scm)
  basis_obj <- BPRMeth::create_rbf_object(M = 4)
  
  set.seed(15)
  # Run Melissa with K = 4 clusters
  melissa_obj <- Melissa::melissa(X = melissa_obj$met, K = 4, basis = basis_obj,
                                  vb_max_iter = 20, vb_init_nstart = 1, 
                                  is_parallel = FALSE)
  
  return(scm)
}

generate_melissa_object <- function (scm, maxgap) {
  
  . <- NULL
  
  # Convert Granges to genomic interval [-1,1]
  chrom_size <- scm@metadata$chrom_sizes
  chrom_start <- sapply(coverage(scm), function(x) {x@lengths[1]}+1)
  interval <- as.data.table(rowRanges(scm))[,.(seqnames,start)]
  
  for (i in 1: length(chrom_size)) {
    interval[seqnames==names(chrom_start)[i], interval := round((2*((start-chrom_start[i])/chrom_size[i]))-1, digits = 3)]
  }
  
  mcols(scm) <- cbind(mcols(scm),interval = interval[,interval])
  
  # Put into Melissa met format
  loc <- paste0(rowRanges(scm)@seqnames,":",chrom_start)
  cells <- list()
  
  for (n in 1:ncol(scm)) {
    
    met <- matrix(cbind(mcols(scm)$interval,get_matrix(scm[,n],assay="binary")),
                  ncol = 2,dimnames =list(loc,NULL))
    
    ends <- len <- seqnames(rowRanges(scm))@lengths
    for (i in 1:length(ends)) ends[i] <- sum(as.vector(len[1:i]))
    starts <- head(c(1, ends + 1), -1)
    
    mets <- lapply(1:length(starts), function (i) as.matrix(met[starts[i]:ends[i],,drop=FALSE]))
    mets <- lapply(mets, function (m) m[-which(m[,2]==-1),])
    
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
  return(melissa_obj)
}

#------------------------------------------------------------------------------------------------------------
#' Imputes the an \code{\link{scMethrix}} object using a iterative PCA model
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object.
#' @param scm A \code{\link{scMethrix}} object
#' @param assay the assay to impute
#' @param new_assay the new assay name
#' @param n_pc The number of principal components to use. This can be a range - e.g. c(1,5) - or a single value. 
#' For a range, the optimal value will be estimated; this is time-intensive.
#' @param ... Additional arguments for missMDA::imputePCA
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' impute_by_iPCA(scMethrix_data, assay = "score", new_assay = "impute", n_pc = 2)
#' @export
#' @references Bro, R., Kjeldahl, K. Smilde, A. K. and Kiers, H. A. L. (2008) Cross-validation of component models: A critical look at current methods. Analytical and Bioanalytical Chemistry, 5, 1241-1251.
#' @references Josse, J. and Husson, F. (2011). Selecting the number of components in PCA using cross-validation approximations. Computational Statistics and Data Analysis. 56 (6), pp. 1869-1879.
impute_by_iPCA <- function(scm = NULL, assay = "score", new_assay = NULL, n_pc = 2, ...) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }
  
  if (new_assay %in% SummarizedExperiment::assayNames(scm)) {
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)
  }
  
  if (is_h5(scm)) {
    stop("HDF5-based objects cannot be imputed by iPCA", call. = FALSE)
  }
  
  if (length(n_pc) == 2) {
    warning("Caution: n_pc is given as range. This can be very time-intensive.")
    n_pc <- missMDA::estim_ncpPCA(get_matrix(scm,assay = assay),ncp.min = n_pc[1], ncp.max = n_pc[2], 
                         method.cv = "Kfold", verbose = TRUE)
    n_pc <- n_pc$ncp
  }
  
  impute <- missMDA::imputePCA(get_matrix(scm,assay = assay), ncp = n_pc, ...)
  assays(scm)[[new_assay]] <- as(impute$completeObs,class(get_matrix(scm,assay=assay))[[1]])
  
  return(scm)
}
