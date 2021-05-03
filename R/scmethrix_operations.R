#' Extract and summarize methylation or coverage info by regions of interest
#' @details Takes \code{\link{scMethrix}} object and summarizes regions
#' @param m \code{\link{scMethrix}} object
#' @param regions genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GenomicRanges}} object
#' @param type matrix which needs to be summarized. Could be `M`, `C`. Default 'M'
#' @param how mathematical function by which regions should be summarized. Can be one of the following: mean, sum, max, min. Default 'mean'
#' @param overlap_type defines the type of the overlap of the CpG sites with the target region. Default value is `within`. For detailed description,
#' see the \code{findOverlaps} function of the \code{\link{IRanges}} package.
#' #param elementMetadata.col columns in \code{\link{scMethrix}}@elementMetadata which needs to be summarised. Default = NULL.
#' @param n_chunks Number of chunks to split the \code{\link{scMethrix}} object in case it is very large. Default = 1.
#' #param n_cores Number of parallel instances. \code{n_cores} should be less than or equal to \code{n_chunks}. If \code{n_chunks} is not specified, then \code{n_chunks} is initialized to be equal to \code{n_cores}. Default = 1.
#' @param verbose Boolean to output progress messages. Default TRUE
#' @param group a column name from sample annotation that defines groups. In this case, the number of min_samples will be
#' tested group-wise.
#' @return table of summary statistic for the given region
#' @examples
#' data('scMethrix_data')
#' get_region_summary(m = scMethrix_data$mem,
#' regions = data.table(chr = c('chr1','chr2'), start = c(1,5), end =  c(5,10)),
#' type = 'M', how = 'mean')
#' @export
get_region_summary = function (m, regions = NULL, n_chunks=1, n_threads = 1, type="M", how = "mean", overlap_type = "within",
                               verbose = TRUE, group = NULL) {

  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!is.null(group) && !(group %in% colnames(m@colData))){
    stop(paste("The column name ", group, " can't be found in colData. Please provid a valid group column."))
  }
  
  if (n_chunks > nrow(m)) {
    n_chucks <- nrow(m)
    warning("n_chunks exceeds number of files. Defaulting to n_chunks = ",n_chunks)
  }
  
  if(verbose) message("Generating region summary...",start_time())
  
  type = match.arg(arg = type, choices = c('M', 'C'))
  how = match.arg(arg = how, choices = c('mean', 'median', 'max', 'min', 'sum', 'sd'))
  yid  <- NULL
  
  regions = cast_granges(regions)
  regions$rid <- paste0("rid_", 1:length(regions))

  overlap_indices <- as.data.table(GenomicRanges::findOverlaps(m, regions, type = overlap_type)) #GenomicRanges::findOverlaps(rowRanges(m), regions)@from
  
  if(nrow(overlap_indices) == 0){
    stop("No overlaps detected")
  }
  
  if(nrow(overlap_indices) == 0){
    stop("No overlaps detected")
  }
  
  colnames(overlap_indices) <- c("xid", "yid")
  overlap_indices[,yid := paste0("rid_", yid)]
  n_overlap_cpgs = overlap_indices[, .N, yid]
  colnames(n_overlap_cpgs) = c('rid', 'n_overlap_CpGs')
  
  if(n_chunks==1){
    if (type == "M") {
      dat = get_matrix(m = m[overlap_indices$xid,], type = "M", add_loci = TRUE)
    } else if (type == "C") {
      dat = get_matrix(m = m[overlap_indices$xid,], type = "C", add_loci = TRUE)
    }
  } else {
    
    stop("Chunking not enabled")
    
    # if(nrow(overlap_indices) < n_chunks){
    #   n_chunks <- nrow(overlap_indices)
    #   warning("Fewer overlaps indicies than n_chunks. Defaulting to n_chunks = ",n_chunks)
    # }
    # 
    # if (n_chunks < n_threads) {
    #   n_threads <- n_chunks
    #   warning("n_threads < n_chunks. Defaulting to n_threads = ",n_threads)
    # }
    # 
    # cl <- parallel::makeCluster(n_threads)  
    # doParallel::registerDoParallel(cl)  
    # 
    # parallel::clusterEvalQ(cl, c(library(data.table), sink(paste0("D:/Git/scMethrix/", Sys.getpid(), ".txt"))))
    # parallel::clusterExport(cl,list('m','scMethrix','type', 'get_matrix','start_time','split_time','stop_time'))
    # 
    # chunk_overlaps <- split(overlap_indices$xid, ceiling(seq_along(overlap_indices$xid) /
    #                                                     ceiling(length(overlap_indices$xid)/n_chunks)))
    # 
    # data = c(parallel::parLapply(cl,chunk_overlaps,fun=function(i) {
    #   cat("Looking for i =",i,typeof(i))
    #   get_matrix(m[i,], type = type, add_loci = TRUE) # TODO: object of type 'S4' is not subsettable
    # }))
    # 
    # # data = lapply(chunk_overlaps,FUN=function(i) {
    # #   cat("Looking for i =",toString(i))
    # #   get_matrix(m[i,], type = type, add_loci = TRUE)
    # # })
    # 
    # parallel::stopCluster(cl)
    # 
    # data <- rbindlist(data)
     
  }
  
  if(nrow(overlap_indices) != nrow(dat)){
    stop("Something went wrong")
  }
  
  dat = cbind(overlap_indices, dat)
  
  #message("-Summarizing overlaps..\n")
  if(how == "mean") {
    message("Summarizing by average")

    output = dat[, lapply(.SD, mean, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(m)))]
  } else if (how == "median") {
    message("Summarizing by median")
    output = dat[, lapply(.SD, median, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(m)))]
  } else if (how == "max") {
    message("Summarizing by maximum")
    output = suppressWarnings(
      dat[, lapply(.SD, max, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(m)))])

    for (j in 1:ncol(output)) set(output, which(is.infinite(output[[j]])), j, NA)
  } else if (how == "min") {
    message("Summarizing by minimum")
    output = suppressWarnings(

      dat[, lapply(.SD, min, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(m)))])
    for (j in 1:ncol(output)) set(output, which(is.infinite(output[[j]])), j, NA)
  } else if (how == "sum") {
    message("Summarizing by sum")
    output = dat[, lapply(.SD, sum, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(m)))]
  }
  
  output = merge(regions, output, by.x = 'rid', by.y = 'yid', all.x = TRUE)

  output = merge(n_overlap_cpgs, output, by = 'rid')
  output$rid <- as.numeric(gsub("rid_","",output$rid))
  
  output <- output[order(output$rid),]
  setnames(output, "seqnames", "chr")
  
  keep <- c("chr", "start", "end", "n_overlap_CpGs", "rid", colnames(m))
  output <- output[, keep, with=FALSE]
  
  if(verbose) message("Region summary generating in ",stop_time())
  
  return(output)
}


#--------------------------------------------------------------------------------------------------------------------------

#' Filter matrices by coverage
#' @details Takes \code{\link{scMethrix}} object and filters CpGs based on coverage statistics
#' @param m \code{\link{scMethrix}} object
#' @param cov_thr minimum coverage required to call a loci covered
#' @param min_samples Minimum number of samples that should have a loci with coverage >= \code{cov_thr}. If \code{group} is given, then this applies per group. Only need one of \code{prop_samples} or \code{min_samples}.
#' @param prop_samples Minimum proportion of samples that should have a loci with coverage >= \code{cov_thr}. If \code{group} is given, then this applies per group. Only need one of \code{prop_samples} or \code{min_samples}.
#' @param group a column name from sample annotation that defines groups. In this case, the number of min_samples will be
#' tested group-wise.
#' @param n_chunks Number of chunks to split the \code{\link{scMethrix}} object in case it is very large. Default = 1.
#' @param n_threads Number of parallel instances. \code{n_threads} should be less than or equal to \code{n_chunks}. If \code{n_threads} is not specified, then \code{n_chunks} is initialized to be equal to \code{n_threads}. Default = 1.
#' @importFrom methods is as new
#' @examples
#' @return An object of class \code{\link{scMethrix}}
#' @export
coverage_filter <- function(m, cov_thr = 1, min_samples = NULL, prop_samples=NULL, group = NULL, n_chunks=1, n_threads=1) {
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!is.null(min_samples) && !is.null(prop_samples)) {
    warning("Both min_samples and prop_samples set. Defaulting to min_samples")
    prop_samples <- NULL
  } else if (is.null(min_samples) && is.null(prop_samples)) {
    stop("Neither min_samples and prop_samples is defined")
  }
  
  if (!is.numeric(cov_thr)) stop("cov_thr is not numeric") 
  
  if (!(is.numeric(min_samples) | is.numeric(prop_samples))){
    stop("min_samples and prop_samples variables must be numeric or NULL")
  }
  
  if (!is.null(group) && !(group %in% colnames(m@colData))){
    stop(paste("The column name ", group, " can't be found in colData. Please provid a valid group column."))
  }
  
  if (n_threads > n_chunks){
    n_chunks <- n_threads
    message("n_threads should be set to be less than or equal to n_chunks.", "\n", "n_chunks has been set 
            to be equal to n_threads = ", n_threads)
  }
  
  message("Filtering by coverage...",start_time())
  
  if (is_h5(m)) {
    if (n_chunks == 1) {
      cov_dat = get_matrix(m = m, type = "C")
      if (!is.null(group)) {
        stop("Groups not implemented yet")
        # row_idx <- sapply(unique(m@colData[, group]), function(c) {
        #   res <- DelayedMatrixStats::rowSums2(cov_dat[, m@colData[, 
        #                                                           group] == c] >= cov_thr, na.rm = TRUE)
        #   row_idx <- (res >= max(min_samples, ceiling(prop_samples * 
        #                                                 sum(m@colData[, group] == c))))
        # })
        # row_idx <- DelayedMatrixStats::rowAlls(row_idx)
      } else {
        res <- DelayedMatrixStats::rowSums2(cov_dat >= cov_thr, na.rm = TRUE)
        row_idx <- (res >= max(min_samples, ceiling(prop_samples * 
                                                      ncol(cov_dat))))
        # row_idx <- (res >= (if (is.null(min_samples)) min_samples else ceiling(prop_samples * 
        #  
      }
    } else {
      message("Only 1 chunk supported for now")
      # 
      # 
      # if (!is.null(group)) {
      #   row_idx <- unlist(mclapply(mc.cores = n_cores, 1:n_chunks, 
      #                              function(i) {
      #                                cov_dat = get_matrix(m[((i - 1) * ceiling(nrow(m)/n_chunks) + 1):min(i * ceiling(nrow(m)/n_chunks), nrow(m)), ], 
      #                                                     type = "C")
      #                                row_idx <- sapply(unique(m@colData[, group]), function(c) {
      #                                  res <- DelayedMatrixStats::rowSums2(cov_dat[, m@colData[,group] == c] >= cov_thr, na.rm = TRUE)
      #                                  row_idx <- (res >= max(min_samples, ceiling(prop_samples * sum(m@colData[, group] == c))))
      #                                })
      #                                row_idx <- DelayedMatrixStats::rowAlls(row_idx)
      #                              }))
      # } else {
      #   row_idx <- unlist(mclapply(mc.cores = n_cores, 1:n_chunks, 
      #                              function(i) {
      #                                cov_dat = get_matrix(m[((i - 1) * ceiling(nrow(m)/n_chunks) + 1):min(i * ceiling(nrow(m)/n_chunks), nrow(m)), ], 
      #                                                     type = "C")
      #                                res <- DelayedMatrixStats::rowSums2(cov_dat >= cov_thr, na.rm = TRUE)
      #                                row_idx <- (res >= max(min_samples, ceiling(prop_samples * ncol(cov_dat))))
      #                              }))
      # }
    }
  } else {
    cov_dat = get_matrix(m = m, type = "C")
    if (!is.null(group)) {
      stop("Groups not implemented yet")
      # row_idx <- sapply(unique(m@colData[, group]), function(c) {
      #   res <- matrixStats::rowSums2(cov_dat[, m@colData[, group] == 
      #                                          c] >= cov_thr, na.rm = TRUE)
      #   row_idx <- (res >= max(min_samples, ceiling(prop_samples * sum(m@colData[, group] == c))))
      # })
      # row_idx <- matrixStats::rowAlls(row_idx)
    } else {
      res <- matrixStats::rowSums2(cov_dat >= cov_thr, na.rm = T)
      row_idx <- (res >= max(min_samples, ceiling(prop_samples * ncol(cov_dat))))
    }
  }
  
  gc()
  message(paste0("Retained ", format(sum(row_idx), big.mark = ","),
                 " of ", format(nrow(m), big.mark = ","), " sites"))
  
  message("Filtered in ",stop_time())
  
  return(m[row_idx, ])
  
}

#--------------------------------------------------------------------------------------------------------------------------
#' Saves an HDF5 \code{\link{scMethrix}} object
#' @details Takes \code{\link{scMethrix}} object and saves it in the specified directory
#' @param m \code{\link{scMethrix}} object
#' @param h5_dir The directory to use. Created, if not existing. Default NULL
#' @param replace Should it overwrite the pre-existing data? FALSE by default.
#' @param ... Parameters to pass to saveHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @examples
#' data('scMethrix_data')
#' dir <- paste0(tempdir(),"/h5")
#' m <- convert_methrix(scMethrix_data$mem, h5_dir=dir)
#' save_HDF5_scMethrix(m, h5_dir = dir, replace = TRUE)
#' @return Nothing
#' @export
save_HDF5_scMethrix <- function(m, h5_dir = NULL, replace = FALSE, ...) {
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is.null(dir)) {
    stop("Please provide the target directory containing ")
  }
  
  if (is_h5(m)) {
    HDF5Array::saveHDF5SummarizedExperiment(x = m, dir = h5_dir, replace = replace, ...)
  } else {
    stop("The object is not a methrix object or not in an HDF5 format. ")
  }
}

#--------------------------------------------------------------------------------------------------------------------------
#' Loads HDF5 scMethrix object
#' @details Takes  directory with a previously saved HDF5Array format \code{\link{scMethrix}} object and loads it
#' @param dir The directory to read in from. Default NULL
#' @param ... Parameters to pass to \code{\link{loadHDF5SummarizedExperiment}}
#' @return An object of class \code{\link{scMethrix}}
#' @examples
#' data('scMethrix_data')
#' dir <- paste0(tempdir(),"/h5")
#' m <- convert_methrix(scMethrix_data$mem, h5_dir=dir)
#' save_HDF5_scMethrix(m, h5_dir = dir, replace = TRUE)
#' n <- load_HDF5_scMethrix(dir)
#' @export
load_HDF5_scMethrix <- function(dir = NULL, ...) {
  
  if (is.null(dir)) {
    stop("Please provide the target directory containing ")
  }
  
  message("Loading HDF5 object", start_time())
  
  m <- HDF5Array::loadHDF5SummarizedExperiment(dir = dir, ...)
  m <- as(m, "scMethrix")
  
  message("Loaded in ",stop_time())
  
  return(m)
}


#--------------------------------------------------------------------------------------------------------------------------
#' Converts HDF5 \code{\link{scMethrix}} object to standard in-memory object.
#' @details Takes a \code{\link{scMethrix}} object and returns with the same object with in-memory assay slots.
#' @param m An object of class \code{\link{scMethrix}}, HDF5 format
#' @return An object of class \code{\link{scMethrix}}
#' @importFrom SummarizedExperiment assays assays<-
#' @examples
#' data('scMethrix_data')
#' dir <- paste0(tempdir(),"/h5")
#' m <- convert_methrix(scMethrix_data$mem, h5_dir=dir)
#' convert_HDF5_methrix(m)
#' @export
convert_HDF5_methrix <- function(m = NULL) {
  
  if (is.null(m) | !is(m, "scMethrix")) {
    stop("Input must be of type scMethrix.")
  }
  if (!is_h5(m)) {
    stop("Input scMethrix must be in HDF5 format.")
  }
  
  message("Converting in-memory scMethrix to HDF5")#,start_time())
  
  assays(m)[[1]] <- as.matrix(assays(m)[[1]])
  m@metadata$is_h5 <- FALSE
  
  # message("Converted in ",stop_time())
  
  return(m)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Converts an in-memory object to an on-disk HDF5 object.
#' @details Takes a \code{\link{scMethrix}} object and returns with the same object with delayed array assay slots
#' with HDF5 backend. Might take long time!
#' @param m An object of class \code{\link{scMethrix}}
#' @param h5_dir Directory for the HDF5 object to be stored
#' @param verbose flag to output messages or not
#' @return An object of class \code{\link{scMethrix}}, HDF5 format
#' @importFrom SummarizedExperiment assays
#' @examples
#' data('scMethrix_data')
#' convert_methrix(scMethrix_data$mem, h5_dir=paste0(tempdir(),"/h5"))
#' @export
convert_methrix <- function(m = NULL, h5_dir = NULL, verbose = TRUE) {
  
  if (is.null(m) | !is(m, "scMethrix")) {
    stop("Input must be of type scMethrix.")
  }
  if (is_h5(m)) {
    stop("Input scMethrix is already in HDF5 format.")
  }
  
  if (verbose) message("Converting in-memory scMethrix to HDF5")#, start_time())
  
  m <- create_scMethrix(methyl_mat = assays(m)[[1]], h5_dir = h5_dir,
                        rowRanges = rowRanges(m), is_hdf5 = TRUE, genome_name = m@metadata$genome,
                        colData = m@colData, chrom_sizes = m@metadata$chrom_sizes, 
                        desc = m@metadata$descriptive_stats, replace = TRUE, verbose = verbose)
  
  #if (verbose) message("Converted in ", stop_time())
  
  return(m)
}



#--------------------------------------------------------------------------------------------------------------------------
#' Order scMethrix object by SD
#' @details Takes \code{\link{scMethrix}} object and reorganizes the data by standard deviation
#' @param m \code{\link{scMethrix}} object
#' @param zero.rm Removes zero values from equations (the default empty value for sparse matrices)
#' @param na.rm Removes the NA values from equations
#' @return An object of class \code{\link{scMethrix}}
#' @examples
#' @export
order_by_sd <- function (m, zero.rm = FALSE, na.rm = FALSE) {
  # 
  # if (!is(m, "scMethrix")){
  #   stop("A valid scMethrix object needs to be supplied.")
  # }
  # 
  # sds <- DelayedMatrixStats::rowSds(x = get_matrix(m),na.rm=TRUE)
  # 
  # m$sd <- sds
  
}

#--------------------------------------------------------------------------------------------------------------------------

#' Subsets an \code{\link{scMethrix}} object based on given conditions.
#' @details Takes \code{\link{scMethrix}} object and filters CpGs based on region, contig and/or sample. Can 
#' either subset to or filter out the input parameters.
#' @param m \code{\link{scMethrix}} object
#' @param regions genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param contigs string of chromosome names to subset by
#' @param samples string of sample names to subset by
#' @param by string to decide whether to "include" or "exclude" the given criteria from the subset
#' @param verbose flag to output messages or not
#' @importFrom IRanges subsetByOverlaps
#' @examples
#' data('scMethrix_data')
#' 
#' contigs <- c("chr1","chr3")
#' regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,5)) 
#' samples <- c("df1","df3")

#' #Subset to only samples bed1 and bed3, and chromosome 1
#' subset_scMethrix(scMethrix_data$mem, samples = samples, contigs = contigs, by = "include")
#' 
#' #Subset to only region "chr1:1-5"
#' subset_scMethrix(scMethrix_data$mem, regions = regions, by = "include")
#' 
#' #Subset to exclude samples bed1 and bed3, and chromosome 1
#' subset_scMethrix(scMethrix_data$mem, samples = samples, contigs = contigs, by = "exclude")
#' 
#' #Subset to exclude region "chr1:1-5"
#' subset_scMethrix(scMethrix_data$mem, regions = regions, by = "exclude")
#' @return An object of class \code{\link{scMethrix}}
#' @export

subset_scMethrix <- function(m, regions = NULL, contigs = NULL, samples = NULL, by=c("include","exclude"), verbose=TRUE) {
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is.null(regions) & is.null(contigs) & is.null(samples)) {
    warning("At least 1 argument mandatory for subsetting. No subset generated")
    return(m)
  }
  
  if (by == "exclude") {
    
    if (!is.null(regions)) {
      regions <- cast_granges(regions)
      reg <- subsetByOverlaps(rowRanges(m), regions, invert = TRUE, type="any", maxgap=-1L, minoverlap=0L)
      m <- subset_scMethrix(m,regions=reduce(reg),by="include")
    }
    
    if (!is.null(contigs)) {
      c <- as.character(seqnames(m)@values)
      m <- subset_scMethrix(m,contigs = c[!c %in% contigs],by="include")
    }
    
    if (!is.null(samples)) {
      s <- row.names(colData(m))
      m <- subset_scMethrix(m,samples = s[!s %in% samples],by="include")
    }
    
  } else {
    
    if (verbose) message("Subsetting CpG sites...",start_time())
    
    if (!is.null(regions)) {
      regions <- cast_granges(regions)
      if (verbose) message("   Subsetting by regions")
      subset <- cast_granges(regions)
      m <- m[GenomicRanges::findOverlaps(rowRanges(m), regions)@from]
      if (nrow(m) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
    }
    
    if (!is.null(contigs)) {
      if (verbose) message("   Subsetting by contigs")
      m <- subset(m, subset = as.vector(seqnames(m)) %in% contigs)
      if (nrow(m) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
    }
    
    if (!is.null(samples)) {
      if (verbose) message("   Subsetting by samples")
      m <- subset(m, select = row.names(colData(m)) %in% samples)
      if (length(m) == 0) stop("Samples not present in the object", call. = FALSE)
    }
    
    if (verbose) message("Subset in ",stop_time())
  }
  
  return(m)
  
}

#--------------------------------------------------------------------------------------------------------------------------
#' Estimate descriptive statistics
#' @details Calculate descriptive statistics
#' @param m \code{\link{scMethrix}} object
#' @param per_chr Estimate stats per chromosome. Default TRUE
#' @examples
#' data('scMethrix_data')
#' 
#' #Get stats for each sample and chromosome
#' get_stats(scMethrix_data$mem)
#' 
#' #Get stats for each sample
#' get_stats(scMethrix_data$mem,per_chr = FALSE)
#' @return data.table of summary stats
#' @export

get_stats <- function(m, per_chr = TRUE) {
  
  if (!is(m, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  message("Getting descriptive statistics...",start_time())
  
  ends <- seqnames(m)@lengths
  for (i in 1:length(ends))
    ends[i] <- sum(as.vector(ends[1:i]))
  starts <- head(c(1, ends + 1), -1)
  
  if (is_h5(m)) {
    if (per_chr) {
      stats <- lapply(1:length(starts), function(x) {
        data.table::data.table(
          Chr = levels(seqnames(m))[x],
          Sample_Name = colnames(m),
          mean_meth = DelayedMatrixStats::colMeans2(get_matrix(m), rows = starts[x]:ends[x], na.rm = TRUE),
          median_meth = DelayedMatrixStats::colMedians(get_matrix(m), rows = starts[x]:ends[x], na.rm = TRUE),
          sd_meth = DelayedMatrixStats::colSds(get_matrix(m), rows = starts[x]:ends[x], na.rm = TRUE)
        )
      })
      
      stats <- data.table::rbindlist(l = stats, use.names = TRUE)
      
    } else {
      stats <- data.table::data.table(
        Sample_Name = colnames(m),
        mean_meth = DelayedMatrixStats::colMeans2(get_matrix(m), na.rm = TRUE),
        median_meth = DelayedMatrixStats::colMedians(get_matrix(m), na.rm = TRUE),
        sd_meth = DelayedMatrixStats::colSds(get_matrix(m), na.rm = TRUE)
      )
    }
    
  } else {
    if (per_chr) {
      stats <- lapply(1:length(starts), function(x) {
        data.table::data.table(
          Chr = levels(seqnames(m))[x],
          Sample_Name = colnames(m),
          mean_meth = matrixStats::colMeans2(get_matrix(m),rows = starts[x]:ends[x],na.rm = TRUE),
          median_meth = matrixStats::colMedians(get_matrix(m), rows = starts[x]:ends[x], na.rm = TRUE),
          sd_meth = matrixStats::colSds(get_matrix(m),rows = starts[x]:ends[x],na.rm = TRUE)
        )
      })
      stats <- data.table::rbindlist(l = stats, use.names = TRUE)
      
    } else {
      stats <-
        data.table::data.table(
          Sample_Name = colnames(m),
          mean_meth = matrixStats::colMeans2(get_matrix(m), na.rm = TRUE),
          median_meth = matrixStats::colMedians(get_matrix(m), na.rm = TRUE),
          sd_meth = matrixStats::colSds(get_matrix(m), na.rm = TRUE)
        )
    }
  }
  
  gc()
  message("Finished in ", stop_time())
  
  return(stats)
}

#--------------------------------------------------------------------------------------------------------------------------

#' Extract methylation or coverage matrices
#' @details Takes \code{\link{scMethrix}} object and returns the \code{methylation} matrix
#' @param m \code{\link{scMethrix}} object
#' @param add_loci Default FALSE. If TRUE adds CpG position info to the matrix and returns as a data.table
#' @param in_granges Do you want the outcome in \code{\link{GRanges}}?
#' @param type Which matrix to get, "m": methlation, "c": coverage
#' @return Coverage or Methylation matrix
#' @examples
#' data('scMethrix_data')
#' 
#' # Get methylation data
#' get_matrix(scMethrix_data$mem)
#' 
#' # Get methylation data with loci
#' get_matrix(scMethrix_data$mem, add_loci=TRUE)
#' 
#' # Get methylation data with loci inside a Granges object 
#' get_matrix(scMethrix_data$mem, add_loci=TRUE, in_granges=TRUE)
#' @export
get_matrix <- function(m, add_loci = FALSE, in_granges=FALSE, type = "M") {
  
  if (!is(m, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (add_loci == FALSE & in_granges == TRUE) {
    warning("Without genomic locations (add_loci= FALSE), it is not possible to convert the results to GRanges, ", 
            "the output will be a data.frame object. ")
  }
  
  type = match.arg(arg = type, choices = c('M', 'C'))
  
  type <- if (type == "M") 1 else 2

  mtx <- SummarizedExperiment::assay(x = m, i = type)
  
  if (add_loci) {
    
    if (is_h5(m)) mtx <- as.data.frame(mtx)
    
    mtx <- as.data.frame(cbind(as.data.frame(rowRanges(m))[,1:3], mtx))
    
    if (in_granges) {
      mtx <- GenomicRanges::makeGRangesFromDataFrame(mtx, keep.extra.columns = TRUE)
      
    } else {
      data.table::setDT(x = mtx)
      colnames(mtx)[1] <- "chr" #TODO: figure out why seqnames is used instead of chr
    }
    
  }
  
  return (mtx)
}

#' Remove loci that are uncovered across all samples
#' @details Takes \code{\link{scMethrix}} object and removes loci that are uncovered across all samples
#' @param m \code{\link{scMethrix}} object
#' @return An object of class \code{\link{scMethrix}}
#' @examples
#' data('scMethrix_data')
#' # Remove uncovered CpGs after subsetting to a single sample
#' remove_uncovered(subset_scMethrix(scMethrix_data$mem, samples = "df1", by="include"))
#' @export
remove_uncovered <- function(m) {
  
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  message("Removing uncovered CpGs...", start_time())
  
  row_idx <- rowSums(!is.na(get_matrix(m)))==0
  
  message(paste0("Removed ", format(sum(row_idx), big.mark = ","),
                 " [", round(sum(row_idx)/nrow(m) * 100, digits = 2), "%] uncovered loci of ",
                 format(nrow(m), big.mark = ","), " sites (",stop_time(),")"))
  
  if (!sum(row_idx) == 0) m <- m[!row_idx, ]
  
  return(m)
}



#--------------------------------------------------------------------------------------------------------------------------

#' Masks too high or too low counts
#' @details Takes \code{\link{scMethrix}} object and masks sites with too high or too low coverage
#'  by putting NA for coverage and beta value. The sites will remain in the object. 
#' @param m \code{\link{scMethrix}} object
#' @param low_count The minimal coverage allowed. Everything below, will get masked. Default = NULL, nothing gets masked.
#' @param high_quantile The quantile limit of coverage. Quantiles are calculated for each sample and everything that belongs to a
#' higher quantile than the defined will be masked. Default = 0.99.
#' @param n_threads Number of parallel instances. Can only be used if \code{\link{scMethrix}} is in HDF5 format. Default = 1.
#' @param type Whether to use the "coverage" matrix or sample "count" when masking
#' @return An object of class \code{\link{scMethrix}}
#' @importFrom SummarizedExperiment assays assays<-
#' @examples
#' @export

mask_methrix <- function(m, low_count = 0, high_quantile = NULL, n_threads=1 ,type="count") {
  
  if (!is(m, "scMethrix")) stop("A valid scMethrix object needs to be supplied.")
  
  if (!is_h5(m) & n_threads != 1) 
     stop("Parallel processing not supported for a non-HDF5 scMethrix object due to probable high memory usage. \nNumber of cores (n_threads) needs to be 1.")
   
  type = match.arg(arg = type, choices = c('count', 'coverage'))
  
  if (!has_cov(m) && type == "coverage") stop("No coverage matrix is present in the object. 
                                              Retry with type='count'")
  
  if (!is.null(high_quantile)) {
    if (type == 'count') stop("high_quantile cannot be used with 'count' ")
    if (high_quantile >= 1 | high_quantile <= 0) stop("High quantile should be between 0 and 1. ")
  }
  
  message("Masking CpG sites in ",high_quantile," quintile and ",type, " < ",low_count, start_time())
  
  if (!is.null(low_count)) {
    
    if(!is.numeric(low_count)){
      stop("low_count must be a numeric value.")
    }
  
    n <- nrow(m) - DelayedMatrixStats::colCounts(get_matrix(m), value = as.integer(NA))
    
    if (type == "count") {
      row_idx <- DelayedMatrixStats::rowCounts(get_matrix(m), 
                                               value = as.integer(NA)) > (nrow(colData(m))-low_count)
    } else {
      row_idx <- DelayedMatrixStats::rowSums2(get_matrix(m,type="C"),na.rm=TRUE) < low_count
    }
    
    if (sum(row_idx) == 0) stop("No CpGs found with low_count")
    
    row_idx <- which(row_idx)  
    emp <- array(as.integer(NA),c(length(row_idx), nrow(colData(m))))

    assays(m)[[1]][row_idx,] <- emp
    if (type == "coverage") assays(m)[[2]][row_idx,] <- emp 
    
    n <- n-(nrow(m) - DelayedMatrixStats::colCounts(get_matrix(m), value = as.integer(NA)))
    
    for (i in seq_along(colnames(m))) {
      if (n[i]==0) next
      message(paste0("   Masked ", n[i], " CpGs due to too low count in sample ",  colnames(m)[i], "."))
    }
  }

  if (!is.null(high_quantile)) {

    message("Masking coverage higher than ",high_quantile * 100," percentile")

    quantiles <- DelayedMatrixStats::colQuantiles(assays(m)[[2]], probs = high_quantile,
                                               na.rm = TRUE, drop = F)
    quantiles <- as.vector(quantiles)
    names(quantiles) <- rownames(m@colData)
    
    row_idx2 <- t(t(assays(m)[[2]]) > quantiles)
    
    if (sum(row_idx2, na.rm = TRUE) == 0) stop("No samples found within in the quantile")
    
    assays(m)[[1]][row_idx2] <- as.integer(NA)
    assays(m)[[2]][row_idx2] <- as.integer(NA)
    
    n <- DelayedMatrixStats::colSums2(row_idx2, na.rm = T)
    
    for (i in seq_along(colnames(m))) {
      if (n[i]==0) next
      message(paste0("Masked ", n[i], " CpGs due to too high coverage in sample ",
                     colnames(row_idx2)[i], "."))
    }
  }
  
  message("Masked in ",stop_time())
  
  return(m)
}
