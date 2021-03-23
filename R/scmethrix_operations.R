get_region_summary = function (m) {
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
}



#--------------------------------------------------------------------------------------------------------------------------
#' Loads HDF5 scMethrix object
#' @details Takes  directory with a previously saved HDF5Array format \code{\link{scMethrix}} object and loads it
#' @param dir The directory to read in from. Default NULL
#' @param ... Parameters to pass to loadHDF5SummarizedExperiment
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data('scMethrix_data')
#' methrix_data_h5 <- convert_methrix(m=methrix_data)
#' target_dir = paste0(getwd(), '/temp1/')
#' save_HDF5_methrix(methrix_data_h5, dir = target_dir, replace = TRUE)
#' load_HDF5_methrix(target_dir)
#' @export
load_HDF5_scMethrix <- function(dir = NULL, ...) {
  
  if (is.null(dir)) {
    stop("Please provide the target directory containing ")
  }
  
  m <- HDF5Array::loadHDF5SummarizedExperiment(dir = dir, ...)
  m <- as(m, "scMethrix")
  return(m)
}


#--------------------------------------------------------------------------------------------------------------------------
#' Converts HDF5 methrix object to standard in-memory object.
#' @details Takes a \code{\link{methrix}} object and returns with the same object with in-memory assay slots.
#' @param m An object of class \code{\link{methrix}}, HDF5 format
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data(methrix_data)
#' m2 <- convert_methrix(m=methrix_data)
#' m <- convert_HDF5_methrix(m=m2)
#' @export
convert_HDF5_methrix <- function(m = NULL) {
  
  if (is.null(m) | !is(m, "scMethrix")) {
    stop("Input must be of type scMethrix.")
  }
  if (!is_h5(m)) {
    stop("Input scMethrix must be in HDF5 format.")
  }

  assays(m)[[1]] <- as.matrix(assays(m)[[1]])
  m@metadata$is_h5 <- FALSE
  return(m)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Converts an in-memory object to an on-disk HDF5 object.
#' @details Takes a \code{\link{methrix}} object and returns with the same object with delayed array assay slots
#' with HDF5 backend. Might take long time!
#' @param m An object of class \code{\link{methrix}}
#' @return An object of class \code{\link{methrix}}, HDF5 format
#' @examples
#' data(methrix_data)
#' m2 <- convert_methrix(m=methrix_data)
#' @export
convert_methrix <- function(m = NULL, h5_dir = NULL) {
  
  if (is.null(m) | !is(m, "scMethrix")) {
    stop("Input must be of type scMethrix.")
  }
  if (is_h5(m)) {
    stop("Input scMethrix is already in HDF5 format.")
  }
  
  m <- create_scMethrix(methyl_mat = assays(m)[[1]], h5_dir = h5_dir,
                      rowRanges = rowRanges(m), is_hdf5 = TRUE, genome_name = m@metadata$genome,
                      colData = m@colData, chrom_sizes = m@metadata$chrom_sizes, 
                      desc = m@metadata$descriptive_stats)
  return(m)
}


  
#--------------------------------------------------------------------------------------------------------------------------
#' Order methrix object by SD
#' @details Takes \code{\link{scMethrix}} object and reorganizes the data by standard deviation
#' @param m \code{\link{scMethrix}} object
#' @param zero.rm Removes zero values from equations (the default empty value for sparse matrices)
#' @param na.rm Removes the NA values from equations
#' @return An object of class \code{\link{scMethrix}}
#' @examples
#' data('scMethrix_data')
#' order_by_sd(m = scMethrix_data)
#' @export
order_by_sd <- function (m, zero.rm = FALSE, na.rm = FALSE) {
  
#  if (!is(m, "scMethrix")){
#    stop("A valid scMethrix object needs to be supplied.")
 # }

  if (zero.rm) {
    sds <- numeric(nrow(m))
    for (i in 1:length(sds)) {
      sds[i] <- sd(as.numeric(m[i,][m[i,]!=0]), na.rm)
      print(paste(i,": ", sds[i]))
      }
  } else {
    sds <- sparseMatrixStats::rowSds(m,na.rm)
  }
  
  row_order <- order(sds, na.last = TRUE, decreasing = TRUE)
  m <- m[row_order, ]
  
  return (m)
  
}

#--------------------------------------------------------------------------------------------------------------------------

#' Subsets \code{\link{scMethrix}} object based on given conditions.
#' @details Takes \code{\link{scMethrix}} object and filters CpGs based on coverage statistics
#' @param m \code{\link{scMethrix}} object
#' @param regions genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param contigs string of chromosome names to subset by
#' @param samples string of sample names to subset by
#' @examples
#' data('scMethrix_data')
#' contigs <- c("chr1","chr3")
#' regions <- GRanges(seqnames = "chr1", ranges = IRanges(1,5)) 
#' samples <- c("df1","df3")
#' 
#' #Subset to chromosome 1
#' subset_scMethrix(scMethrix_data, contigs = contigs)
#' 
#' #Subset to samples bed1 and bed3
#' subset_scMethrix(scMethrix_data, samples = samples)
#' 
#' #Subset to samples bed1 and bed3, and chromosome 1
#' subset_scMethrix(scMethrix_data, samples = samples, contigs = contigs)
#' 
#' #Subset to region "chr1:1-150"
#' subset_scMethrix(scMethrix_data, regions = regions)

#' @return An object of class \code{\link{scMethrix}}
#' @export

subset_scMethrix <- function(m, regions = NULL, contigs = NULL, samples = NULL, verbose=TRUE) {
  
  message("Subsetting scMethrix")
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is.null(regions) & is.null(contigs) & is.null(samples)) {
    warning("At least 1 argument mandatory for subsetting. No subset generated")
    return(m)
  }
    
  if (!is.null(regions)) {
    message("   Subsetting by regions")
    subset <- cast_granges(regions)
    m <- m[findOverlaps(rowRanges(m), regions)@from]
    if (nrow(m) == 0)
      stop("Subsetting resulted in zero entries", call. = FALSE)
  }
  
  if (!is.null(contigs)) {
    message("   Subsetting by contigs")
    m <- m[seqnames(m) %in% contigs]
    if (nrow(m) == 0)
      stop("Subsetting resulted in zero entries", call. = FALSE)
  }
  
  if (!is.null(samples)) {
    message("   Subsetting by samples")
    m <- subset(m, select = colData(m)[, 1] %in% samples)
    if (length(m) == 0)
      stop("Samples not present in the object", call. = FALSE)
  }

  return(m)
  
}


#--------------------------------------------------------------------------------------------------------------------------
#' Estimate descriptive statistics
#' @details Calculate descriptive statistics
#' @param m \code{\link{methrix}} object
#' @param per_chr Estimate stats per chromosome. Default TRUE
#' @seealso \code{\link{plot_stats}}
#' @examples
#' data('methrix_data')
#' get_stats(methrix_data)
#' @return data.table of summary stats
#' @export
get_stats <- function(m, per_chr = TRUE) {
  median <- . <- sd <- chr <- NULL
  start_proc_time <- proc.time()
  
  if (!is(m, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  #Counting NA's seems to be not-ideal (especially in cases of masked/coverage filtered cases.)
  # row_idx <- data.table::as.data.table(which(is.na(get_matrix(m = m,
  #     type = "C")), arr.ind = TRUE))
  # colnames(row_idx) <- c("row", "col")
  # row_idx <- split(row_idx, as.factor(as.character(row_idx$col)))
  
  if (per_chr) {
    cov_stat <- lapply(1:ncol(m), function(i) {
      get_matrix(m = m[, i],
                 type = "C",
                 add_loci = TRUE)[, c(1, 4), with = FALSE][, .(
                   mean_cov = lapply(.SD,
                                     matrixStats::mean2, na.rm = TRUE),
                   median_cov = lapply(.SD,
                                       median, na.rm = TRUE),
                   sd_cov = lapply(.SD, sd, na.rm = TRUE)
                 ),
                 by = chr]
    })
    
    meth_stat <- lapply(1:ncol(m), function(i) {
      get_matrix(m = m[, i],
                 type = "M",
                 add_loci = TRUE)[, c(1, 4), with = FALSE][, .(
                   mean_meth = lapply(.SD,
                                      matrixStats::mean2, na.rm = TRUE),
                   median_meth = lapply(.SD,
                                        median, na.rm = TRUE),
                   sd_meth = lapply(.SD, sd, na.rm = TRUE)
                 ),
                 by = chr]
    })
    
    names(meth_stat) <-  colnames(m)
    names(cov_stat) <- colnames(m)
    
    cov_stat <-
      data.table::rbindlist(l = cov_stat,
                            use.names = TRUE,
                            idcol = "Sample_Name")
    meth_stat <-
      data.table::rbindlist(l = meth_stat,
                            use.names = TRUE,
                            idcol = "Sample_Name")
    stats <-
      merge(meth_stat, cov_stat, by = c("chr", "Sample_Name"))
    colnames(stats)[1] <- "Chromosome"
    stats$Chromosome <-
      factor(x = stats$Chromosome,
             levels = m@metadata$chrom_sizes$contig)
  } else {
    if (is_h5(m)) {
      stats <- data.table::data.table(
        Sample_Name = colnames(m),
        mean_cov = DelayedMatrixStats::colMeans2(get_matrix(m = m, "C"), na.rm = TRUE),
        median_cov = DelayedMatrixStats::colMedians(get_matrix(m = m, "C"), na.rm = TRUE),
        sd_cov = DelayedMatrixStats::colSds(get_matrix(m = m, "C"), na.rm = TRUE)
      )
      
    } else {
      stats <-
        data.table::data.table(
          Sample_Name = colnames(m),
          mean_meth = matrixStats::colMeans2(get_matrix(m = m, "M"), na.rm = TRUE),
          median_meth = matrixStats::colMedians(get_matrix(m = m, "M"), na.rm = TRUE),
          sd_meth = matrixStats::colSds(get_matrix(m = m, "M"), na.rm = TRUE)
        )
    }
  }
  
  gc()
  message("-Finished in:  ", data.table::timetaken(start_proc_time))
  
  return(stats)
}


#--------------------------------------------------------------------------------------------------------------------------

#' Extract methylation or coverage matrices
#' @details Takes \code{\link{methrix}} object and returns user specified \code{methylation} or \code{coverage} matrix
#' @param m \code{\link{methrix}} object
#' @param type can be \code{M} or \code{C}. Default 'M'
#' @param add_loci Default FALSE. If TRUE adds CpG position info to the matrix and returns as a data.table
#' @param in_granges Do you want the outcome in \code{GRanges}?
#' @return Coverage or Methylation matrix
#' @examples
#' data('methrix_data')
#' #Get methylation matrix
#' get_matrix(m = methrix_data, type = 'M')
#' #Get methylation matrix along with loci
#' get_matrix(m = methrix_data, type = 'M', add_loci = TRUE)
#' #' #Get methylation data as a GRanges object
#' get_matrix(m = methrix_data, type = 'M', add_loci = TRUE, in_granges=TRUE)
#' @export
get_matrix <- function(scm, add_loci = FALSE, in_granges=FALSE) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (add_loci == FALSE & in_granges == TRUE) {
    warning("Without genomic locations (add_loci= FALSE), it is not possible to convert the results to GRanges, ", 
            "the output will be a data.frame object. ")
  }

  mtx <- as.data.frame(SummarizedExperiment::assay(x = scm, i = 1))
  
  if (add_loci) {
    
    mtx <- as.data.frame(cbind(as.data.frame(rowRanges(scm))[,1:3],
                               mtx))
      
    if (in_granges) mtx <- GenomicRanges::makeGRangesFromDataFrame(mtx, keep.extra.columns = TRUE)

  }

  return (mtx)
}
  
  
  
  
  
  
  
  
