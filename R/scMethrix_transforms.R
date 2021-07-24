#------------------------------------------------------------------------------------------------------------
#' Transforms an assay in an \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object. It is typically used for The transform is applied column-wise to optimize how HDF5 files access sample data. If HDF5 objects are used, transform functions should be from \pkg{DelayedMatrixStats}.
#' @inheritParams generic_scMethrix_function
#' @param h5_temp string; temporary directory to store the temporary HDF5 files
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
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object. Typically, most assays will use either mean (for measurements) or sum (for counts). The transform is applied column-wise to optimize how HDF5 files access sample data. If HDF5 objects are used, transform functions should be  from \pkg{DelayedMatrixStats}.
#' 
#' In the output object, the number of CpGs in each region is saved in mcol(scm)$n_cpgs.
#' @inheritParams generic_scMethrix_function
#' @param regions The regions from which to make the bins
#' @param bin_size integer; The size of each bin. First bin will begin at the start position of the first genomic
#' region on the chromosome. If NULL, there will be one bin per region. Default 100000.
#' @param bin_by character; can create bins by # of base pairs "bp" or by # of CpG sites "cpg". Default "bp"
#' @param trans named vector of closures; The transforms for each assay in a named vector. Default NULL, meaning that 
#' operations for "counts" assay is sum(x, na.rm=TRUE), and for all other assays is mean(x, na.rm=TRUE)
#' @param overlap_type character; defines the type of the overlap of the CpG sites with the target region. 
#' Default value is `within`. For detailed description, see the \code{findOverlaps} function of the 
#' \code{\link{IRanges}} package.
#' @param h5_dir directory to store an H5 based object
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' regions <- GRanges(seqnames = c("chr1"), ranges = IRanges(1,200000000)) 
#' regions <- unlist(tile(regions,10))
#' bin_scMethrix(scMethrix_data, regions = regions)
#' @export
bin_scMethrix <- function(scm = NULL, regions = NULL, bin_size = 100000, bin_by = "bp", trans = NULL, 
                          overlap_type = "within", h5_dir = NULL, verbose = TRUE, n_chunks = 1, n_threads = 1) {

  yid <- NULL
  
  if (is.null(trans)) {
    trans <- c(counts = function(x) sum(x,na.rm=TRUE))
  }
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
 # if (is_h5(scm) && is.null(h5_dir)) stop("Output directory must be specified", call. = FALSE)
  
  bin_by = match.arg(arg = bin_by, choices = c("bp","cpg"))
  
  if (!is.null(regions)) {
    regions = cast_granges(regions)
    scm <- subset_scMethrix(scm, regions = regions) 
  } else { # If no region is specifed, use entire chromosomes
    regions = range(rowRanges(scm))
  }
  
  if (verbose) message("Subsetting for ",length(regions)," input regions...",start_time())
  
  regions$rid <- paste0("rid_", 1:length(regions))
  
  overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, regions, type = overlap_type)) 
  
  if(nrow(overlap_indices) == 0){
    stop("No overlaps detected")
  }
  
  colnames(overlap_indices) <- c("xid", "yid")
  overlap_indices[,yid := paste0("rid_", yid)]
  
  if (!is.null(bin_size)) {
    
    if (verbose) message("Binning by ",bin_by," with size of ",bin_size,"...")
    
    if (bin_by == "cpg") {
      
      rrng = GRanges()
      
      for(rid in 1:length(regions$rid)) {
        
        #if (verbose) message("   Processing region ",rid,"/",length(regions$rid))
        
        idx <- which(overlap_indices[,yid == regions$rid[rid]])
        idx <- split_vector(idx,num = bin_size, by = "size")
        
        for(i in idx) {
          gr <- range(rowRanges(scm[c(i[1],i[length(i)]),]))
          gr$n_cpgs <- length(i)
          rrng <- c(rrng,gr)
        }
      }
      
    } else if (bin_by == "bp") {
      
      rrng <- unlist(tile(regions, width = bin_size)) #TODO: Should switch this to using RLE lookup
      
      idx <- as.data.table(GenomicRanges::findOverlaps(scm, rrng, type = overlap_type))
      
      rrng <- rrng[unique(idx$subjectHits)]
      rrng$n_cpgs <- rle(idx$subjectHits)$lengths
    }
    
  } else { # If no bin_size is specified, use the entire region
    rrng <- regions
  }
  
  if (verbose) message("Generated ",length(rrng)," bins in ",split_time())
  
  assays <- list()
  
  rrng$rid <- paste0("rid_", 1:length(rrng))
  overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, rrng, type = overlap_type))
  colnames(overlap_indices) <- c("xid", "yid")
  overlap_indices[,yid := paste0("rid_", yid)]
  
  # Function to process each bin

  for (name in SummarizedExperiment::assayNames(scm)) {

    if (verbose) message("   Filling bins for the ",name," assay...")
    
    
    if (is.null(trans[[name]])) { # If no named vector is specified, default to mean
      op <- function(x) mean(x,na.rm=TRUE)
    } else {
      op <- trans[[name]]
    }
    

    
    if (is_h5(scm)) {
      message("This isn't done yet")
      
      # rid_list <- split_vector(rrng$rid,num = n_threads)
      # 
      # 
      # cl <- parallel::makeCluster(n_threads)  
      # doParallel::registerDoParallel(cl) 
      # parallel::clusterEvalQ(cl, c(library(data.table)))
      # parallel::clusterExport(cl,list('overlap_indices','bin','mtx','yid'), envir = environment())
      # 
    } else {

      worker <- function (mtx,overlap_indices,op) {

        mtx <- mtx[,lapply(.SD,op),by=(overlap_indices$yid)]
        mtx <- mtx[,overlap_indices:=NULL]
        return(mtx)
      }
      
      mtx <- data.table(get_matrix(scm,assay=name))[overlap_indices$xid,] 
      cols <- split_vector(1:ncol(mtx),n_chunks)

      cl <- parallel::makeCluster(n_threads)
      doParallel::registerDoParallel(cl)
      parallel::clusterEvalQ(cl, c(library(data.table)))
      parallel::clusterExport(cl,list('overlap_indices','worker','op'), envir = environment())
      on.exit(parallel::stopCluster(cl))

     
      mtx <- lapply(cols,function(col) mtx[,..col])
      #mtx <- lapply(mtx,worker,overlap_indices = overlap_indices,op = op)
      mtx <- parLapply(cl,mtx,worker,overlap_indices = overlap_indices,op=op)
      
      
      mtx <- setDT(unlist(mtx, recursive = FALSE))

      # Split by rids
      # worker <- function (mtx, yid, op) {
      #   mtx <- mtx[,lapply(.SD,op),by=(yid)]
      #   mtx <- mtx[,yid:=NULL]
      #   return(mtx)
      # }
      # mtx <- data.table(get_matrix(scm,assay=name))[overlap_indices$xid,] 
      # 
      # #Make two chunked lists for xid and yids
      # yid_list <- rle(overlap_indices$yid)
      # yid_list <- split(overlap_indices$yid,ceiling(rep(1:length(yid_list$lengths), yid_list$lengths) /
      #                                                 (length(yid_list$lengths)/n_chunks)))
      # lens <- sapply(yid_list,function(yids) length(yids))      
      # xid_list <- split(overlap_indices$xid,rep(1:length(lens),lens))
      # 
      # cl <- parallel::makeCluster(n_threads)  
      # doParallel::registerDoParallel(cl) 
      # parallel::clusterEvalQ(cl, c(library(data.table)))
      # parallel::clusterExport(cl,list('worker'), envir = environment())
      # on.exit(parallel::stopCluster(cl))
      # 
      # mtx <- clusterMap(cl,worker,mtx=lapply(xid_list,function(xid) mtx[unlist(xid),]), yid = yid_list, MoreArgs = list(op = op))
      # mtx <- rbindlist(lapply(mtx, as.data.frame.list))

      # Mapply algorithm 
      # mtx <- mapply(worker,mtx=lapply(xid_list,function(xid) mtx[unlist(xid),]), yid = yid_list, MoreArgs = list(op = op),simplify=TRUE)
      # cols <- row.names(mtx)
      # mtx <- lapply(1:nrow(mtx),function(x) unlist(mtx[x,]))
      # mtx <- t(rbindlist(lapply(mtx, as.data.frame.list)))
      # colnames(mtx) <- cols
      
      # Basic algorithm
       # mtx <- data.table(get_matrix(scm,assay=name))[overlap_indices$xid,] #TODO: Somehow missing rows if not subset, not sure why
       # mtx <- mtx[,lapply(.SD,op),by=overlap_indices$yid]
       # mtx <- mtx[,overlap_indices:=NULL]
      # 
      assays[[name]] <- as(mtx,class(get_matrix(scm,assay=name)))
      
    }
    
    if (verbose) message("Bins filled in ",split_time())
  }
  
  rrng$rid <- NULL
  
  if (verbose) message("Rebuilding experiment...")
  
  if (is_h5(scm)) {
    m_obj <- create_scMethrix(assays = assays, rowRanges=rrng, is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),
                              replace = replace)  
  } else {
    m_obj <- create_scMethrix(assays = assays, rowRanges=rrng, is_hdf5 = FALSE, 
                              genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),)
  }
  
  if (verbose) message("Experiment binned for ",length(regions)," regions containing ",length(m_obj)," total bins in ",stop_time())
  
  # parallel::stopCluster(cl)
  # rm(cl)
  
  return (m_obj)
}



bin_scMethrix2 <- function(scm = NULL, regions = NULL, bin_size = 100000, bin_by = "bp", trans = NULL, 
                          overlap_type = "within", h5_dir = NULL, verbose = TRUE, n_chunks = 1, n_threads = 1) {
  
  yid <- NULL
  
  if (is.null(trans)) {
    trans <- c(counts = function(x) sum(x,na.rm=TRUE))
  }
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  # if (is_h5(scm) && is.null(h5_dir)) stop("Output directory must be specified", call. = FALSE)
  
  bin_by = match.arg(arg = bin_by, choices = c("bp","cpg"))
  
  if (!is.null(regions)) {
    regions = cast_granges(regions)
    scm <- subset_scMethrix(scm, regions = regions) 
  } else { # If no region is specifed, use entire chromosomes
    regions = range(rowRanges(scm))
  }
  
  if (verbose) message("Subsetting for ",length(regions)," input regions...",start_time())
  
  regions$rid <- paste0("rid_", 1:length(regions))
  
  overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, regions, type = overlap_type)) 
  
  if(nrow(overlap_indices) == 0){
    stop("No overlaps detected")
  }
  
  colnames(overlap_indices) <- c("xid", "yid")
  overlap_indices[,yid := paste0("rid_", yid)]
  
  if (!is.null(bin_size)) {
    
    if (verbose) message("Binning by ",bin_by," with size of ",bin_size,"...")
    
    if (bin_by == "cpg") {
      
      rrng = GRanges()
      
      for(rid in 1:length(regions$rid)) {
        
        #if (verbose) message("   Processing region ",rid,"/",length(regions$rid))
        
        idx <- which(overlap_indices[,yid == regions$rid[rid]])
        idx <- split_vector(idx,num = bin_size, by = "size")
        
        for(i in idx) {
          gr <- range(rowRanges(scm[c(i[1],i[length(i)]),]))
          gr$n_cpgs <- length(i)
          rrng <- c(rrng,gr)
        }
      }
      
    } else if (bin_by == "bp") {
      
      rrng <- unlist(tile(regions, width = bin_size)) #TODO: Should switch this to using RLE lookup
      
      idx <- as.data.table(GenomicRanges::findOverlaps(scm, rrng, type = overlap_type))
      
      rrng <- rrng[unique(idx$subjectHits)]
      rrng$n_cpgs <- rle(idx$subjectHits)$lengths
    }
    
  } else { # If no bin_size is specified, use the entire region
    rrng <- regions
  }
  
  if (verbose) message("Generated ",length(rrng)," bins in ",split_time())
  
  assays <- list()
  
  rrng$rid <- paste0("rid_", 1:length(rrng))
  overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, rrng, type = overlap_type))
  colnames(overlap_indices) <- c("xid", "yid")
  overlap_indices[,yid := paste0("rid_", yid)]
  
  # Function to process each bin
    for (name in SummarizedExperiment::assayNames(scm)) {
      
      if (verbose) message("   Filling bins for the ",name," assay...")
      
      
      if (is.null(trans[[name]])) { # If no named vector is specified, default to mean
        op <- function(x) mean(x,na.rm=TRUE)
      } else {
        op <- trans[[name]]
      }
      
      if (is_h5(scm)) {
        
        rid_list <- split_vector(rrng$rid,num = n_threads)
        
        
        
        
        
        
        
        cl <- parallel::makeCluster(n_threads)  
        doParallel::registerDoParallel(cl) 
        parallel::clusterEvalQ(cl, c(library(data.table)))
        parallel::clusterExport(cl,list('overlap_indices','bin','mtx','yid'), envir = environment())
        
      } else {
        
        mtx <- data.table(get_matrix(scm,assay=name))[overlap_indices$xid,] #TODO: Somehow missing rows if not subset, not sure why
        mtx <- mtx[,lapply(.SD,op),by=overlap_indices$yid]
        mtx <- mtx[,overlap_indices:=NULL]

        assays[[name]] <- as(mtx,class(get_matrix(scm,assay=name)))
        
      }
      
      if (verbose) message("Bins filled in ",split_time())
    }
  
  rrng$rid <- NULL
  
  if (verbose) message("Rebuilding experiment...")
  
  if (is_h5(scm)) {
    m_obj <- create_scMethrix(assays = assays, rowRanges=rrng, is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),
                              replace = replace)  
  } else {
    m_obj <- create_scMethrix(assays = assays, rowRanges=rrng, is_hdf5 = FALSE, 
                              genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),)
  }
  
  if (verbose) message("Experiment binned for ",length(regions)," regions containing ",length(m_obj)," total bins in ",stop_time())
  
  # parallel::stopCluster(cl)
  # rm(cl)
  
  return (m_obj)
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
#' @references Kapourani CA, Sanguinetti G (2019). “Melissa: Bayesian clustering and imputation of single cell methylomes.” Genome Biology, 20, 61. doi: 10.1186/s13059-019-1665-8.
impute_by_melissa <- function (scm, threshold = 50, assay = "score", new_assay = "impute") {
  
  . <- NULL
  
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
    warning("Imputation cannot be done on HDF5 data. Data will be cast as matrix for imputation.")
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
    
    met <- matrix(c(mcols(scm)$interval,as.matrix(get_matrix(scm[,n],assay="binary"))),
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
#' @param ... Additional arguments for \link[missMDA]{imputePCA}
#' @inheritParams impute_regions
#' @inheritParams missMDA::imputePCA
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' impute_by_iPCA(scMethrix_data, assay = "score", new_assay = "impute", n_pc = 2)
#' @export
#' @references Bro, R., Kjeldahl, K. Smilde, A. K. and Kiers, H. A. L. (2008) Cross-validation of component models: A critical look at current methods. Analytical and Bioanalytical Chemistry, 5, 1241-1251.
#' @references Josse, J. and Husson, F. (2011). Selecting the number of components in PCA using cross-validation approximations. Computational Statistics and Data Analysis. 56 (6), pp. 1869-1879.
impute_by_iPCA <- function(scm = NULL, regions = NULL, assay = "score", new_assay = "iPCA", 
                           n_chunks = 1, n_threads = 1, overlap_type = "within", n_pc = 2, ...) {

  if (length(n_pc) > 1) {
    warning("Caution: n_pc is given as range. This can be very time-intensive.")
    n_pc <- missMDA::estim_ncpPCA(as.matrix(get_matrix(scm,assay = assay)),ncp.min = n_pc[1], ncp.max = n_pc[2], 
                                  method.cv = "Kfold", verbose = TRUE)
    n_pc <- n_pc$ncp
  }

  op <- function(mtx) missMDA::imputePCA(mtx, ncp = n_pc, ...)$completeObs
    
  scm <- impute_regions(scm = scm, assay = assay, new_assay = new_assay, regions = regions, op = op, 
                       n_chunks = n_chunks, n_threads = n_threads, overlap_type = overlap_type)
    
  return(scm)
    
}

#------------------------------------------------------------------------------------------------------------
#' Imputes the an \code{\link{scMethrix}}object using a random forest model
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object.
#' @param ... Additional arguments for \link[missForest]{missForest}
#' @inheritParams impute_regions
#' @inheritParams missForest::missForest
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' impute_by_RF(scMethrix_data, assay = "score", new_assay = "impute")
#' @export
#' @import missForest
#' @references Stekhoven, D. J., & Bühlmann, P. (2012). MissForest—non-parametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.
impute_by_RF <- function(scm = NULL, regions = NULL, assay = "score", new_assay = "RF", 
                         n_chunks = 1, n_threads = 1, overlap_type = "within", ...) {

  # impute <- missForest::missForest(as.matrix(get_matrix(scm,assay = assay)))#, ...)
  # assays(scm)[[new_assay]] <- as(impute$ximp,class(get_matrix(scm,assay=assay))[[1]])
  
  op <- function(mtx) missForest::missForest(mtx, ...)$ximp
  
  scm <- impute_regions(scm = scm, assay = assay, new_assay = new_assay, regions = regions, op = op, 
                     n_chunks = n_chunks, n_threads = n_threads, overlap_type = overlap_type)

  return(scm)
}


#------------------------------------------------------------------------------------------------------------
#' Imputes the an \code{\link{scMethrix}}object using a random forest model
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object.
#' @param ... Additional arguments for \link[impute]{impute.knn}
#' @inheritParams impute_regions
#' @inheritParams impute::impute.knn
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' impute_by_kNN(scMethrix_data)
#' @export
#' @import impute
#' @references Hastie T, Tibshirani R, Narasimhan B, Chu G (2021). impute: impute: Imputation for microarray data. R package version 1.66.0.
impute_by_kNN <- function(scm = NULL, regions = NULL, assay = "score", new_assay = "kNN", 
                          n_chunks = 1, n_threads = 1, overlap_type = "within", k = 10, ...) {
  
  op <- function(mtx) impute::impute.knn(mtx, k = min(k,ncol(mtx)), 
                           rowmax = 1.0, colmax = 1.0, maxp = 1500, ...)$data
  
  scm <- impute_regions(scm = scm, assay = assay, new_assay = new_assay, regions = regions, op = op, 
                     n_chunks = n_chunks, n_threads = n_threads, overlap_type = overlap_type)

  # impute <- impute::impute.knn(as.matrix(get_matrix(scm,assay = assay)), k = min(k,ncol(scm)), 
  #                      rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=123)
  # 
  # assays(scm)[[new_assay]] <- as(impute$data,class(get_matrix(scm,assay=assay))[[1]])
  # 
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Generic imputation return function
#' @details Uses the specified imputation operation to evaluation an scMethrix object.
#' @param regions Granges; the regions to impute. Default is by chromosome.
#' @param overlap_type string; 
#' @param op closure; the imputation operation
#' @inheritParams generic_scMethrix_function
#' @return list; two \code{\link{scMethrix}} objects names 'training' and 'test'
#' @examples
#' data('scMethrix_data')
#' generate_training_set(scMethrix_data, training_prop = 0.2)
#' @export
impute_regions <- function(scm = NULL, assay="score", new_assay = NULL, regions = NULL, op = NULL, n_chunks = 1, 
                               n_threads = 1, overlap_type="within") {
  
  yid <- NULL
  
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
    warning("Imputation cannot be done on HDF5 data. Data will be cast as matrix for imputation.")
  }
  
  if (!is.null(regions)) {
    regions = cast_granges(regions)
    scm <- subset_scMethrix(scm, regions = regions) 
  } else { # If no region is specifed, use entire chromosomes
    regions = range(rowRanges(scm))
  }
  
  assays(scm)[[new_assay]] <- assays(scm)[[assay]]
  
  regions$rid <- paste0("rid_", 1:length(regions))
  overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, regions, type = overlap_type))
  colnames(overlap_indices) <- c("xid", "yid")
  overlap_indices[,yid := paste0("rid_", yid)]
  rid_list <- split_vector(regions$rid,num = n_chunks)
  rid_list <- lapply(rid_list, function(rids) {split_vector(rids,num = n_threads)})

  for (rids in rid_list) {
   
    idx <- lapply(rids,FUN = function(rid) overlap_indices[overlap_indices$yid %in% rid]$xid)
    impute <- lapply(idx,function(i) op(get_matrix(scm[i,],assay))) # This is in parallel later
    
    for (i in 1:length(idx)) assays(scm)[[new_assay]][idx[[i]],] <- impute[[i]]
    
  }
   
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Splits an scMethrix object into two for use as a training and test set
#' @details Typically used for teaching classification algorithms. The seed can be set for consistency.
#' @param training_prop numeric; The size of the training set as a proportion of the experiment (0 to 1)
#' For a range, the optimal value will be estimated; this is time-intensive.
#' @param seed string; value to use for sampling
#' @inheritParams generic_scMethrix_function
#' @return list; two \code{\link{scMethrix}} objects names 'training' and 'test'
#' @examples
#' data('scMethrix_data')
#' generate_training_set(scMethrix_data, training_prop = 0.2)
#' @export
generate_training_set <- function(scm = NULL, training_prop = 0.2, seed = "123") {
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (training_prop > 1 || training_prop < 0) {
    stop("training_prop must in the range of [0,1]", call. = FALSE)
  }
  set.seed(seed)
  idx <- sort(sample(1:nrow(scm),floor(nrow(scm)*training_prop)))
  
  training <- scm[idx,]
  test <- scm[setdiff(1:nrow(scm),idx),]
  
  return(list(training = training,test = test))
}

#------------------------------------------------------------------------------------------------------------
#' Generates a random subset of CpG sites
#' @details From an \code{\link{scMethrix}} object, this will randomly select \code{n_cpgs} and create a new
#' object containing only those CpGs. This is typically used for approximation or visualization. The seed
#' can be specified for consistency. 
#' @param n_cpgs numeric; The number of CpGs to include
#' @param seed string; value to use for sampling
#' @inheritParams generic_scMethrix_function
#' @return scMethrix; an experiment with n_cpgs
#' @examples
#' data('scMethrix_data')
#' generate_random_subset(scMethrix_data,n_cpgs = round(nrow(scMethrix_data)/2))
#' @export
generate_random_subset <- function(scm = NULL, n_cpgs = 10000, seed = "123") {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (n_cpgs > nrow(scm) || n_cpgs < 1) {
    n_cpgs = max(1,n_cpgs)
    n_cpgs = min(nrow(scm),n_cpgs)
    warning("Invalid n_cpgs. Must be between 1 and ",nrow(scm),". Defaulted to ",n_cpgs)
  }
  
  set.seed(seed)
  idx <- sort(sample(1:nrow(scm),n_cpgs))
  
  return(scm[idx,])
}