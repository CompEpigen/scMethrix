
#' Ternary operator
#' @details Acts in the same way as ternary operators in other language. Must surround with brackets to use)
#' @param x The value to return if true
#' @param y The value to return if false
#' @return The x or y value depending on boolean
#' @import
#' @examples
`?` <- function(x, y)
  eval(
    sapply(
      strsplit(
        deparse(substitute(y)), 
        ":"
      ), 
      function(e) parse(text = e)
    )[[2 - as.logical(x)]])

#' Checks if scMethrix object is an HDF5 object
#' @details Acts in the same way as ternary operators in other language. Must surround with brackets to use)
#' @param m The scMethrix object
#' @return boolean Whether the object is HDF5
#' @import
#' @examples
is_h5 = function(m) {
  return(m@metadata$is_h5)
}

#' Returns file name minus the extension from a file.path to represent the sample name
#' @details Acts in the same way as ternary operators in other language. Must surround with brackets to use)
#' @param s A file.path
#' @return string containing the sample name
#' @import tools
#' @examples
#' get_sample_name("c:/folder/folder/sample.name.ext")
get_sample_name = function(s) {
  
  return(tools::file_path_sans_ext(basename(s)))
}

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
chunk_granges = function(gr,factor = NA, percent = NA, num = NA) { #=NULL, percent = NULL
  
  if (length(which(is.na(c(factor,percent,num))))!=2) stop("1 argument mandatory for chunking.")
  
  if (!is.na(num)) {
    splits <- length(gr)%/%num
    splits <- num*(1:splits-1)
  }
  
  if (!is.na(percent)) factor = 100/percent
  
  if (!is.na(factor)) {
    num <- floor(length(gr)/factor)
    splits <- floor(length(gr)/num)
    splits <- 1+rep(0:(splits-1))*num
  }

  grl <- List()
  
  for (i in 1:length(splits)) {
    s <- splits[i]
    grl[[i]] <- gr[s:(s+num-1)]
  }
    
  grl[[i+1]] <- gr[(last(splits)+num):length(gr)]  
  
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
#' @export
start_time <- function() {
  assign("time.all", proc.time()["elapsed"], envir=baseenv())
  assign("time.split", proc.time()["elapsed"], envir=baseenv())
  invisible(NULL)
}

#' Outputs the split/lap/iteration time 
#' @details Gets the stored elapsed proc.time() from either the initial start_time or the previous split_time
#' @return Returns formatted elapsed time since start_time or last split_time
#' @export
split_time <- function() {
  
  time <- get("time.split", envir=baseenv())
  if (is.na(time)) stop("start_time() not set")
  time <- proc.time()["elapsed"]-time
  assign("time.split", proc.time()["elapsed"], envir=baseenv())
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}

#' Stops an internal stopwatch and outputs overall time
#' @details Gets the stored elapsed proc.time() from initial start_time() to calculate overall runtime
#' @return Returns formatted elapsed time since start_time
#' @export
stop_time <- function() {
  
  time <- get("time.all", envir=baseenv())
  if (is.na(time)) stop("start_time() not set")
  time <- proc.time()["elapsed"]-get("time.all", envir=baseenv())
  assign("time.split", NA, envir=baseenv())
  assign("time.all", NA, envir=baseenv())
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}



