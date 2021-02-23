get_region_summary = function (m) {
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
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
#' regions <- GRanges(seqnames = "chr1", ranges = IRanges(1,150)) 
#' samples <- c("bed1","bed3")
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
  
  if (is_ondisk(m)) {
    
    files <- get_files(m)
    subset <- rowRanges(m)
    
    if (length(files) == 0) stop("No tabix files to subset", call. = FALSE)
    
    if (!is.null(samples)) {
      message("   Subsetting by samples")
      idx <- unlist(lapply(files,get_sample_name))
      idx <- which(idx %in% samples)   
      files <- files[idx]
      if (length(files) == 0) stop("Samples not present in files", call. = FALSE)
    }
    
    if (!is.null(contigs)) {
      message("   Subsetting by contigs")
      subset <- subset[seqnames(subset) == contigs]
      if (length(subset) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
    }
    
    if (!is.null(regions)) {
      message("   Subsetting by regions")
      regions <- cast_granges(regions)
      subset <- subsetByOverlaps(subset, regions)
      if (length(subset) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
    }
    
    subset <- sort(subset) # TODO: Somehow subset becomes unsorted; not sure why
    overlaps <- NULL
    
    for (file in files) {
      score <- get_tabix_scores(file, subset)
      if (is.null(overlaps)) {overlaps <- score
      } else {overlaps <- cbind(overlaps,score)}
    }
    
    if (nrow(overlaps) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
    
    m_obj <- overlaps
    
  } else {
  
    if (!is.null(regions)) {
      message("   Subsetting by regions")
      regions <- cast_granges(regions)
      overlaps <- subsetByOverlaps(m, regions) 
      if (nrow(overlaps) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
      m_obj <- overlaps
    }
    
    if (!is.null(contigs)) {
      message("   Subsetting by contigs")
      overlaps <- m[seqnames(m) == contigs]
      if (nrow(overlaps) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
      m_obj <- overlaps
    }
    
    if (!is.null(samples)) {
      message("   Subsetting by samples")
      overlaps <- m
      for (sample in samples) {mcols(overlaps)[[sample]] <- NULL}
      if (length(overlaps) == 0) stop("Samples not present in the object", call. = FALSE)
      m_obj <- overlaps
    }
  }
  
  return(m_obj)
  
}
  