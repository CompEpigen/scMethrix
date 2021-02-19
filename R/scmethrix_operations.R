library("sparseMatrixStats")
library("rbenchmark")


get_region_summary = function (m) {
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
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
  
subset_methrix= function () {}
  
  
  
  
coverage_filter = function () {}


get_matrix = function {}
  
  
scMethrix2bsseq 


remove_uncovered



region_filter


mask_methrix



combine_methrix




#--------------------------------------------------------------------------------------------------------------------------

#' Subsets \code{\link{scMethrix}} object based on given conditions.
#' @details Takes \code{\link{scMethrix}} object and filters CpGs based on coverage statistics
#' @param m \code{\link{scMethrix}} object
#' @param regions genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param contigs chromosome names to subset by
#' @param samples sample names to subset by
#' @examples
#' data('scMethrix_data')
#' #Subset to chromosome 1
#' subset_scMethrix(scMethrix_data, contigs = 'chr1')
#' @return An object of class \code{\link{scMethrix}}
#' @export
subset_scMethrix() <- function(m, regions = NULL, contigs = NULL, samples = NULL) {
  
  if (!is(m, "scMmethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!is.null(regions)) {
  
    target_regions <- cast_ranges(regions)
    overlaps <- subsetByOverlaps(m, target_regions)  
    
    if (nrow(overlaps) == 0) {
      stop("Subsetting resulted in zero entries")
    }
    
    m <- overlaps
  }
  
  if (!is.null(contigs)) {
    
    contigs  <- GRanges(seqnames=contigs)
    overlaps <- subsetByOverlaps(m, contigs)
    
    if (nrow(overlaps) == 0) {
      stop("Subsetting resulted in zero entries")
    }
    
    m <- overlaps
  }
  
  if (!is.null(samples)) {
    message("Subsetting by samples")
    
    overlaps <- m
    
    for (sample in samples) {mcols(overlaps)[[sample]] <- NULL}
      
    if (length(overlaps) == 0) {
      stop("None of the samples are present in the object")
    }
    
    m <- overlaps
  }
  
  return(m)
  
}

srapply <- function(s, func) {
  
  
mcols(data.mg)  
  
  
  mcols(gr)$value <- NULL
  
  
  
}


  