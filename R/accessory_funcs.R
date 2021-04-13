#' Checks if scMethrix object is an HDF5 object
#' @details Acts in the same way as ternary operators in other language. Must surround with brackets to use)
#' @param m The scMethrix object
#' @return boolean Whether the object is HDF5
#' @examples
#' data('scMethrix_data')
#' is_h5(scMethrix_data$mem)
#' @export
is_h5 = function(m) {
  return(m@metadata$is_h5)
}

#' Returns file name minus the extension from a file.path to represent the sample name
#' @details Acts in the same way as ternary operators in other language. Must surround with brackets to use)
#' @param s A file.path
#' @return string containing the sample name
#' @import tools
#' @examples
#' #get_sample_name("C:/dir/dir/filename.ext")
 #' @export
get_sample_name = function(s) {
  
  return(tools::file_path_sans_ext(basename(s)))
}

#' Splits a vector into subvectors
#' @details Splits a vector into n chunks or chunks of n size.
#' @param v The vector to split
#' @param n The number to split by
#' @param by Whether to split by chunk or by size
#' @return vector
#' @examples
#' # Split vector into 4 sub vectors
#' split_vector(c(1,2,3,4,5,6,7,8),4,by="chunk")
#' 
#' # Split vector into sub-vectors with a size of 2
#' split_vector(c(1,2,3,4,5,6,7,8),2,by="size")
#' @export
split_vector = function(v,n, by = c("chunk","size")) {
  
  if (match.arg(by)=="size") {
    
    return (split(v, ceiling(seq_along(v)/n)))
    
  } else {
  
    x <- length(v)
    
    if (x/n < 2) stop("Length of input vector must be at least 2x greater than the number of chunks")
    
    e <- ceiling(x/n)*(1:(n-1))
    s <- c(0,e+1)
    e <- c(e,x)
    
    return(lapply(1:n, function(i) {
      return(v[s[i]:e[i]])
    }))
  }
}

#' Chunks a Granges object by factor, percent or number
#' @details Divides Granges into a list based on a chunking parameter
#' @param gr The Granges object
#' @param factor The factor to which divide the chunks into
#' @param percent The percentage of the Granges to chunk
#' @param num The number of regions to include in each chunk
#' @return GRangesList containing all the chunked Granges
#' @import GenomicRanges
#' @examples
#' @export
chunk_granges = function(gr,factor = NA, percent = NA, num = NA) { #=NULL, percent = NULL
  
  if (sum(is.na(c(factor,percent,num))) != 2) stop("Max 1 argument for chunking.")
  
  if (!is.na(percent)) factor = 100/percent
  
  if (!is.na(factor)) {
    num <- floor(length(gr)/factor)
  }
  
  if (!is.na(num)) {
    splits <- ceiling(length(gr)/num)
    splits <- num*(0:(splits-1))+1
  }
  
  grl <- list()
  
  for (i in 1:length(splits)) {
    s <- splits[i]
    grl[[i]] <- gr[s:(s+num-1)]
  }
    
 # grl[[i+1]] <- gr[(last(splits)+num):length(gr)]  
  
  return (GenomicRanges::GRangesList(grl))
}

#--------------------------------------------------------------------------------------------------------------------------
# Parse genomic regions and convert them to key'd data.table
cast_ranges <- function(regions, set.key = TRUE) {
  chr <- . <- NULL
  if (is(regions, "GRanges")) {
    target_regions <- data.table::as.data.table(x = regions)
    target_regions[, `:=`(seqnames, as.character(seqnames))]
    colnames(target_regions)[seq_len(3)] <- c("chr", "start", "end")
    if (set.key){
      data.table::setDT(x = target_regions, key = c("chr", "start", "end"))}
    target_regions <- target_regions[, .(chr, start, end)]
  } else if (is(regions, "data.frame")) {
    if (all(c("chr", "start", "end") %in% colnames(regions))){
      regions <- regions[,c("chr", "start", "end")]
    } else {
      warning("Columns with names chr, start and end are not found. Assuming that the first three columns are chr, start and end.")
    }
    target_regions <- data.table::as.data.table(x = regions)
    colnames(target_regions)[seq_len(3)] <- c("chr", "start", "end")
    target_regions <- target_regions[, .(chr, start, end)]
    target_regions[, `:=`(chr, as.character(chr))]
    target_regions[, `:=`(start, as.numeric(as.character(start)))]
    target_regions[, `:=`(end, as.numeric(as.character(end)))]
    if (set.key){
      data.table::setDT(x = target_regions, key = c("chr", "start", "end"))}
  } else {
    stop("Invalid input class for regions. Must be a data.table or GRanges object")
  }
  
  target_regions
}

#------------------------------------------------------------------------
# Replace region (text or otherwise) with Granges object
cast_granges <- function(regions) {
  if (is(regions, "GRanges")) return (regions)
  
}

#' Starts an internal stopwatch
#' @details Save the current time to later use for split/lap and overall times
#' @return NULL
start_time <- function() {
  assign("time.all", proc.time()["elapsed"], envir=baseenv())
  assign("time.split", proc.time()["elapsed"], envir=baseenv())
  invisible(NULL)
}

#' Outputs the split/lap/iteration time 
#' @details Gets the stored elapsed proc.time() from either the initial start_time or the previous split_time
#' @return Returns formatted elapsed time since start_time or last split_time
split_time <- function() {
  
  time <- get("time.split", envir=baseenv())
  if (!is.numeric(time)) stop("start_time() not set")
  time <- proc.time()["elapsed"]-time
  assign("time.split", proc.time()["elapsed"], envir=baseenv())
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}

#' Stops an internal stopwatch and outputs overall time
#' @details Gets the stored elapsed proc.time() from initial start_time() to calculate overall runtime
#' @return Returns formatted elapsed time since start_time
stop_time <- function() {
  
  time <- get("time.all", envir=baseenv())
  if (!is.numeric(time)) stop("start_time() not set")
  time <- proc.time()["elapsed"]-get("time.all", envir=baseenv())
  assign("time.split", NA, envir=baseenv())
  assign("time.all", NA, envir=baseenv())
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}

#' Subsets a given list of CpGs by another list of CpGs
#' @details Typically used to reduce the number of potential CpG sites to include only those present 
#' in the input files so as to maximize performance and minimize resources. Can also be used for quality
#' control to see if there is excessive number of CpG sites that are not present in the reference genome.
#' @param ref_cpgs A reference set of CpG sites (e.g. Hg19 or mm10) in bedgraph format
#' @param gen_cpgs A subset of CpG sites. Usually obtained from read_index.
#' @return Returns list of CpG sites in bedgraph format
subset_ref_cpgs <- function(ref_cpgs, gen_cpgs, verbose = TRUE) {
  keys <- plyr::join.keys(ref_cpgs, gen_cpgs, c("chr","start"))
  sub_cpgs <- ref_cpgs[keys$x %in% keys$y, , drop = FALSE]
  r <- nrow(ref_cpgs)
  s <- nrow(sub_cpgs)
  g <- nrow(gen_cpgs)
  if (verbose) message("Dropped ",r-s,"/",r," CpGs (",round((r-s)/r*100,2),"%) from the reference set")
  if (verbose) message(g-s,"/",g," subset CpGs (",round((g-s)/g*100,2),"%) were not present in the reference set")
  return(sub_cpgs)
}

