get_region_summary = function (m) {
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
}

#--------------------------------------------------------------------------------------------------------------------------
#' Saves an HDF5 \code{\link{scMethrix}} object
#' @details Takes \code{\link{methrix}} object and saves it in the specified directory
#' @param m \code{\link{methrix}} object
#' @param h5_dir The directory to use. Created, if not existing. Default NULL
#' @param replace Should it overwrite the pre-existing data? FALSE by default.
#' @param ... Parameters to pass to saveHDF5SummarizedExperiment
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
#' @param ... Parameters to pass to loadHDF5SummarizedExperiment
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
#' Converts HDF5 scMethrix object to standard in-memory object.
#' @details Takes a \code{\link{scMethrix}} object and returns with the same object with in-memory assay slots.
#' @param m An object of class \code{\link{scMethrix}}, HDF5 format
#' @return An object of class \code{\link{scMethrix}}
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
#' @return An object of class \code{\link{scMethrix}}, HDF5 format
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
  
#  if (!is(m, "scMethrix")){
#    stop("A valid scMethrix object needs to be supplied.")
 # }
# 
#   if (zero.rm) {
#     sds <- numeric(nrow(m))
#     for (i in 1:length(sds)) {
#       sds[i] <- sd(as.numeric(m[i,][m[i,]!=0]), na.rm)
#       print(paste(i,": ", sds[i]))
#       }
#   } else {
#     sds <- sparseMatrixStats::rowSds(m,na.rm)
#   }
#   
#   row_order <- order(sds, na.last = TRUE, decreasing = TRUE)
#   m <- m[row_order, ]
#   
#   return (m)
  
}

#--------------------------------------------------------------------------------------------------------------------------

#' Filters an \code{\link{scMethrix}} object based on given conditions.
#' @details Takes \code{\link{scMethrix}} object and filters CpGs based on region, contig and/or sample. Can 
#' either subset to or filter out the input parameters.
#' @param m \code{\link{scMethrix}} object
#' @param regions genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param contigs string of chromosome names to subset by
#' @param samples string of sample names to subset by
#' @param verbose flag to output messages or not
#' @examples
#' data('scMethrix_data')
#' 
#' contigs <- c("chr1","chr3")
#' regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,5)) 
#' samples <- c("df1","df3")
#' 
#' #Subset to only chromosome 1
#' subset_scMethrix(scMethrix_data$mem, contigs = contigs, by = "include")
#' 
#' #Subset to only samples bed1 and bed3
#' subset_scMethrix(scMethrix_data$mem, samples = samples, by = "include")
#' 
#' #Subset to only samples bed1 and bed3, and chromosome 1
#' subset_scMethrix(scMethrix_data$mem, samples = samples, contigs = contigs, by = "include")
#' 
#' #Subset to only region "chr1:1-5"
#' subset_scMethrix(scMethrix_data$mem, regions = regions, by = "include")
#' 
#' #Subset to exclude chromosome 1
#' subset_scMethrix(scMethrix_data$mem, contigs = contigs, by = "exclude")
#' 
#' #Subset to exclude samples bed1 and bed3
#' subset_scMethrix(scMethrix_data$mem, samples = samples, by = "exclude")
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
#' @seealso \code{\link{plot_stats}}
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
#' @details Takes \code{\link{scMethrix}} object and returns user specified \code{methylation} or \code{coverage} matrix
#' @param m \code{\link{scMethrix}} object
#' @param type can be \code{M} or \code{C}. Default 'M'
#' @param add_loci Default FALSE. If TRUE adds CpG position info to the matrix and returns as a data.table
#' @param in_granges Do you want the outcome in \code{GRanges}?
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
get_matrix <- function(m, add_loci = FALSE, in_granges=FALSE) {
  
  if (!is(m, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (add_loci == FALSE & in_granges == TRUE) {
    warning("Without genomic locations (add_loci= FALSE), it is not possible to convert the results to GRanges, ", 
            "the output will be a data.frame object. ")
  }

  mtx <- SummarizedExperiment::assay(x = m, i = 1)
  
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
#' @param n_cores Number of parallel instances. Can only be used if \code{\link{scMethrix}} is in HDF5 format. Default = 1.
#' @return An object of class \code{\link{scMethrix}}
#' @examples
#' @export
mask_methrix <- function(m, low_count = NULL, high_quantile = 0.99, n_cores=1) {

  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!is_h5(m) & n_cores != 1) {
    stop("Parallel processing not supported for a non-HDF5 scMethrix object due to probable high memory usage. \nNumber of cores (n_cores) needs to be 1.")
  }
  
  message("Masking CpG sites ( ",high_quintile," quintile and <= ",low_count, " count", start_time())
  
  if (!is.null(low_count)) {
    
    if(!is.numeric(low_count)){
      stop("low_count must be a numeric value.")
    }
    
    message("Masking count lower than ", low_count)

    if(is_h5(m)) {
      n <- nrow(m) - DelayedMatrixStats::colCounts(get_matrix(m), value = as.integer(NA))
      row_idx <- DelayedMatrixStats::rowCounts(get_matrix(m), value = as.integer(NA)) > low_count
      row_idx <- matrix(rep(row_idx,ncol(m)),ncol = ncol(m))
      row_idx <- DelayedArray(row_idx)
      assays(m)[[1]][row_idx] <- as.integer(NA)
      n <- n-(nrow(m) - DelayedMatrixStats::colCounts(get_matrix(m), value = as.integer(NA)))
    } else {
      n <- nrow(m) - MatrixStats::colSums2(is.na(get_matrix(m)))
      row_idx <- !(MatrixStats::rowSums(is.na(get_matrix(m))) <= low_count)
      assays(m)[[1]][!!row_idx,] <- as.integer(NA)
      n <- n - (nrow(m) - MatrixStats::colSums2(is.na(get_matrix(m))))
    }
      
    for (i in seq_along(colnames(m))) {
      message(paste0("Masked ", n[i], " CpGs due to too low count in sample ",  colnames(m)[i], "."))
    }
    
    
  }
  # 
  # 
  # 
  # if (!is.null(high_quantile)) {
  #   if (high_quantile >= 1 | high_quantile <= 0) {
  #     stop("High quantile should be between 0 and 1. ")
  #   }
  #   
  #   message("\n\n-Masking coverage higher than ",
  #           high_quantile * 100,
  #           " percentile")
  #   
  #   
  #   
  #   if (is_h5(m)) {
  #     if (n_cores == 1) {
  #       quantiles <-
  #         DelayedMatrixStats::colQuantiles(assays(m)[[2]],
  #                                          probs = high_quantile,
  #                                          na.rm = TRUE,
  #                                          drop = F)
  #     }
  #     else {
  #       quantiles <-
  #         simplify2array(mclapply(mc.cores = n_cores, 1:ncol(assays(m)[[2]]),
  #                                 function(i)
  #                                   quantile(assays(m)[[2]][, i], probs = high_quantile, na.rm = TRUE)))
  #     }
  #     quantiles <- as.vector(quantiles)
  #     names(quantiles) <- rownames(m@colData)
  #   } else {
  #     quantiles <-
  #       matrixStats::colQuantiles(assays(m)[[2]], probs = high_quantile, na.rm = TRUE)
  #     quantiles <- as.vector(quantiles)
  #     names(quantiles) <- rownames(m@colData)
  #   }
  #   
  #   
  #   row_idx2 <- t(t((assays(m)[[2]])) > quantiles)
  #   assays(m)[[1]][row_idx2] <- as.double(NA)
  #   assays(m)[[2]][row_idx2] <- as.integer(NA)
  #   
  #   
  #   if (is_h5(m)) {
  #     if (n_cores == 1) {
  #       n <- DelayedMatrixStats::colSums2(row_idx2, na.rm = T)
  #     }
  #     else {
  #       n <- simplify2array(mclapply(mc.cores = n_cores, 1:ncol(row_idx2),
  #                                    function(i)
  #                                      sum(row_idx2[, i], na.rm = T)))
  #     }
  #   } else {
  #     n <- colSums(row_idx2, na.rm = T)
  #   }
  #   
  #   
  #   for (i in seq_along(colnames(m))) {
  #     message(paste0(
  #       "-Masked ",
  #       n[i],
  #       " CpGs due to too high coverage in sample ",
  #       colnames(row_idx2)[i],
  #       "."
  #     ))
  #     
  #   }
  #   
  #   
  
  message("Masked in",stop_time())

  return(m)
}
