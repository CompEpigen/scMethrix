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

#' Filter matrices by region
#' @details Takes \code{\link{scMethrix}} object and filters CpGs based on supplied regions in data.table or GRanges format
#' @param m \code{\link{scMethrix}} object
#' @param regions genomic regions to filter-out. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param type defines the type of the overlap of the CpG sites with the target regions. Default value is `any`. For detailed description,
#' see the \code{foverlaps} function of the \code{\link{data.table}} package.
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data('methrix_data')
#' region_filter(m = methrix_data,
#' regions = data.table(chr = 'chr21', start = 27867971, end =  27868103))
#' @export
region_filter <- function(m, regions=NULL, type = "any") {

  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is.null(regions)) {
    warning("No region specified. No filtering performed")
    return(m)
  }

  regions2 <- subsetByOverlaps(rowRanges(m), regions, invert = TRUE, type, maxgap=-1L, minoverlap=0L)
  return(subset_scMethrix(m,regions=reduce(regions2)))
  
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

#TODO: remove the hardlink to the assay
get_stats <- function(m, per_chr = TRUE) {
  median <- . <- sd <- chr <- NULL
  start_proc_time <- proc.time()
  
  if (!is(m, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  ends <- seqnames(m)@lengths
  for (i in 1:length(ends))
    ends[i] <- sum(as.vector(ends[1:i]))
  starts <- head(c(1, ends + 1), -1)
  
  if (is_h5(m)) {
    if (per_chr) {
      stats <- lapply(1:length(starts), function(x) {
        data.table::data.table(
          Chr = levels(l@values)[x],
          Sample_Name = colnames(m),
          mean_meth = DelayedMatrixStats::colMeans2(assays(m)[[1]], rows = starts[x]:ends[x], na.rm = TRUE),
          median_meth = DelayedMatrixStats::colMedians(assays(m)[[1]], rows = starts[x]:ends[x], na.rm = TRUE),
          sd_meth = DelayedMatrixStats::colSds(assays(m)[[1]], rows = starts[x]:ends[x], na.rm = TRUE)
        )
      })
      
      stats <- data.table::rbindlist(l = stats, use.names = TRUE)
      
    } else {
      stats <- data.table::data.table(
        Sample_Name = colnames(m),
        mean_meth = DelayedMatrixStats::colMeans2(assays(m)[[1]], na.rm = TRUE),
        median_meth = DelayedMatrixStats::colMedians(assays(m)[[1]], na.rm = TRUE),
        sd_meth = DelayedMatrixStats::colSds(assays(m)[[1]], na.rm = TRUE)
      )
    }
    
  } else {
    if (per_chr) {
      stats <- lapply(1:length(starts), function(x) {
        data.table::data.table(
          Chr = levels(l@values)[x],
          Sample_Name = colnames(m),
          mean_meth = matrixStats::colMeans2(assays(m)[[1]],rows = starts[x]:ends[x],na.rm = TRUE),
          median_meth = matrixStats::colMedians(assays(m)[[1]], rows = starts[x]:ends[x], na.rm = TRUE),
          sd_meth = matrixStats::colSds(assays(m)[[1]],rows = starts[x]:ends[x],na.rm = TRUE)
        )
      })
      stats <- data.table::rbindlist(l = stats, use.names = TRUE)
      
    } else {
      stats <-
        data.table::data.table(
          Sample_Name = colnames(m),
          mean_meth = matrixStats::colMeans2(assays(m)[[1]], na.rm = TRUE),
          median_meth = matrixStats::colMedians(assays(m)[[1]], na.rm = TRUE),
          sd_meth = matrixStats::colSds(assays(m)[[1]], na.rm = TRUE)
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
get_matrix <- function(m, add_loci = FALSE, in_granges=FALSE) {
  
  if (!is(m, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (add_loci == FALSE & in_granges == TRUE) {
    warning("Without genomic locations (add_loci= FALSE), it is not possible to convert the results to GRanges, ", 
            "the output will be a data.frame object. ")
  }

  mtx <- as.data.frame(SummarizedExperiment::assay(x = m, i = 1))
  
  if (add_loci) {
    
    mtx <- as.data.frame(cbind(as.data.frame(rowRanges(m))[,1:3], (is_h5(m) ? as.data.frame(mtx) : mtx)))
      
    if (in_granges) {
      mtx <- GenomicRanges::makeGRangesFromDataFrame(mtx, keep.extra.columns = TRUE)

    } else {
      data.table::setDT(x = mtx)
      colnames(mtx)[1] <- "chr" #TODO: figure out why seqnames is used instead of chr
    }
      
  }

  
  return (mtx)
}
  
  
  
  
  
  
  
  
