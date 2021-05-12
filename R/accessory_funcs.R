#' Checks if \code{\link{scMethrix}} object is an HDF5 object
#' @param m The \code{\link{scMethrix}} object
#' @return boolean Whether the object is HDF5
#' @examples
#' data('scMethrix_data')
#' is_h5(scMethrix_data$mem)
#' @export
is_h5 = function(m) {
  return(m@metadata$is_h5)
}

#' Checks if \code{\link{scMethrix}} object has a coverage matrix
#' @param m The \code{\link{scMethrix}} object
#' @return boolean Whether the object has a coverage matrix
#' @examples
#' data('scMethrix_data')
#' has_cov(scMethrix_data$mem)
#' @export
has_cov = function(m) {
  return(m@metadata$has_cov)
}

#' Returns file name minus the extension from a file.path to represent the sample name
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
#' @details Splits a vector into \code{num} vectors or vectors of \code{n} size. If \code{len(vec)%%num != 0},
#' then the last vector will have a length of less than size \code{num}
#' @param vec The vector to split
#' @param num The number to split by
#' @param by Whether to split by '\code{chunks}' or by '\code{size}'
#' @return A vector of vectors
#' @examples
#' # Split vector into 4 sub vectors
#' split_vector(c(1,2,3,4,5,6,7,8),4,by="chunk")
#' 
#' # Split vector into sub-vectors with a size of 2
#' split_vector(c(1,2,3,4,5,6,7,8),2,by="size")
#' @export
split_vector = function(vec, num = 1, by = "chunks") {
  
  if (!is.numeric(num)) {
    stop("num must be numeric")
  }
  
  type = match.arg(arg = by, choices = c('chunks', 'size'))
  
  if (type=="size") {
    
    vec <- split(vec, ceiling(seq_along(vec)/num))
    return (unname(vec))
    
  } else {
    
    len <- length(vec)
    
    if (len/num < 2) stop("Length of input vector must be at least 2x greater than the number of chunks")
    
    chunks <- ceiling(len/num)*(1:(num-1))
    idx <- c(0,chunks+1)
    chunks <- c(chunks,len)
    
    return(lapply(1:num, function(i) {
      return(vec[idx[i]:chunks[i]])
    }))
  }
}

#' Chunks a \code{\link{GRanges}} object by \code{factor}, \code{percent} or \code{number}
#' @details Divides \code{\link{GRanges}} into a list based on a chunking parameter
#' @param gr The \code{\link{GRanges}} object
#' @param factor The factor to which divide the chunks into
#' @param percent The size in percentage of CpGs in each chunk
#' @param num The number of CpGs to include in each chunk
#' @return \code{\link{GRangesList}} containing all the chunked \code{\link{GRanges}}
#' @import GenomicRanges
#' @examples
#' @export
split_granges = function(gr,factor = NA, percent = NA, num = NA) { #=NULL, percent = NULL
  
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

#--- cast_granges -------------------------------------------------------------------------------------------
#' Casts regions into \code{\link{GRanges}} format
#' @details Casts the input as a \code{\link{GRanges}} object. Input can be \code{\link{GRanges}} or a 
#' \code{\link{data.frame}}-compatible class that can be cast through \code{as.data.frame()}. Input BED format
#'  must be \code{chr-start-end} for \code{\link{data.frame}} objects.
#' @param regions The input regions
#' @return \code{\link{GRanges}} object with the input regions
#' @import GenomicRanges
#' @examples
#' regions = data.table(chr = 'chr1', start = 1, end = 100)
#' cast_granges(regions) 
#' @export
cast_granges <- function(regions) {
  if (is(regions, "GRanges")) {return (regions)
  } else if (is(regions,"data.frame")) {return (GenomicRanges::makeGRangesFromDataFrame(regions))
  } else {stop("Invalid input class for regions. Must be a GRanges or data.frame-like")}
  return(regions)
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
#' @details Gets the stored elapsed \\code{\link{proc.time}} from either the initial 
#' \code{\link{start_time}} or the previous \code{split_time}
#' @return Returns formatted elapsed time since \code{\link{start_time}} or last \code{\link{split_time}}
#' @export
split_time <- function() {
  time <- get("time.split", envir=baseenv())
  if (!is.numeric(time)) stop("start_time() not set")
  time <- proc.time()["elapsed"]-time
  assign("time.split", proc.time()["elapsed"], envir=baseenv())
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}

#' Stops an internal stopwatch and outputs overall time
#' @details Gets the stored elapsed \code{proc.time()} from initial \code{\link{start_time}} to calculate 
#' overall runtime
#' @return Returns formatted elapsed time since \code{\link{start_time}}
#' @export
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
#' @param gen_cpgs A subset of CpG sites. Usually obtained from \code{\link{read_index}}.
#' @param verbose flag to output messages or not
#' @return Returns list of CpG sites in bedgraph format
#' @export
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

