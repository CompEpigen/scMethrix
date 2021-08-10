#' Checks if \code{\link{scMethrix}} object is an HDF5 object
#' @param m The \code{\link{scMethrix}} object
#' @return boolean Whether the object is HDF5
#' @examples
#' data('scMethrix_data')
#' is_h5(scMethrix_data)
#' @export
is_h5 = function(m) {
  return(m@metadata$is_h5)
}

#' Checks if \code{\link{scMethrix}} object has a coverage matrix
#' @param m The \code{\link{scMethrix}} object
#' @return boolean Whether the object has a coverage matrix
#' @examples
#' data('scMethrix_data')
#' has_cov(scMethrix_data)
#' @export
has_cov = function(m) {
  return("counts" %in% SummarizedExperiment::assayNames(m))
}

#' Returns sample name derived from the input file name
#' @param s A file.path
#' @return string containing the sample name
#' @import tools
#' @examples
#' #get_sample_name("C:/dir/dir/filename.ext")
#' @export
get_sample_name = function(s) {
  if (!is.character(s)) stop("Must be a string file path")
  return(tools::file_path_sans_ext(basename(s)))
}

#' Binarize an input value based on a \code{threshold}
#' @details Assigns a value of 0 or 1 based on being < or > the \code{thresdhold}, respectively.
#'  If \code{x} = \code{threshold}, \code{x} = 0. NA values are assigned as \code{rep_na}
#' @param x A value to binarize
#' @param threshold The threshold for binarizing. Will default to the center number between max and min.
#' @param rep_na The value to replace missing values with. Default NA. 
#' @param verbose boolean; flag for whether to display threshold or not
#' @return 1 or 0, if above of below the threshold, or 'rep.na' if NA
#' @examples
#' vals <- c(0,0.25,0.5,0.75,1,NA)
#' binarize(vals, threshold=0.5)
#' @export
binarize = function(x,threshold = NULL, rep_na = NA, verbose = FALSE) {
  if(is.null(threshold)) {
    threshold <- (max(x,na.rm=TRUE)-min(x,na.rm=TRUE))/2
    if (verbose) message("Threshold: ", threshold)
  }
  x[x <= threshold] <- 0
  x[x  > threshold] <- 1
  x[is.na(x)] <- rep_na
  return(x)
}

#' A faster version of cbind when trying to combine lists of data.tables
#' @details Assigns a value of 0 or 1 based on being < or > the \code{thresdhold}, respectively.
#'  If \code{x} = \code{threshold}, \code{x} = 0. NA values are assigned as \code{rep_na}
#' @param ... A list of data.tables with identical # of rows
#' @return data.table; the cbinded output 
#' @export
colbind = function(...) {
  setDT(
    unlist(..., recursive = FALSE),
    check.names = TRUE
  )[]
}

#' Splits a vector into list of vectors by \code{chunk} or \code{size}
#' @details Splits a vector into a list of \code{num} vectors or vectors of \code{n} size. If \code{len(vec)%%num != 0},
#' then the last vector will have a length of less than size \code{num}
#' @param vec vector; The vector to split
#' @param num integer; The number to split by
#' @param by string; Whether to split by '\code{chunks}' or by '\code{size}'
#' @return A list of vectors
#' @examples
#' # Split vector into 4 sub vectors
#' split_vector(c(1,2,3,4,5,6,7,8),4,by="chunk")
#' 
#' # Split vector into sub-vectors with a size of 2
#' split_vector(c(1,2,3,4,5,6,7,8),2,by="size")
#' @export
split_vector = function(vec, num = 1, by = c('chunks', 'size')) {

  if (!is.numeric(num)) {
    stop("num must be numeric")
  }
  
  type = match.arg(arg = by, choices = c('chunks', 'size'))
  
  if (type=="size") {
    
    vec <- split(vec, ceiling(seq_along(vec)/num))
    return (unname(vec))
    
  } else {
    
    len <- length(vec)

    if (len < num) num <- len
    
    #if (len/num < 2) stop("Length of input vector must be at least 2x greater than the number of chunks")
    
    chunks <- ceiling(len/num)*(1:(num-1))
    idx <- c(0,chunks+1)
    chunks <- c(chunks,len)
    
    return(lapply(1:num, function(i) {
      return(vec[idx[i]:chunks[i]])
    }))
  }
}

#' Splits a \code{\link{GRanges}} object by \code{chunks}, \code{percent} or \code{size}
#' @details Divides each region \code{\link{GRanges}} into equally-sized lists based on a chunking parameter. 
#' If \code{len(gr) %% number != 0}, then the last vector will have a length of less than size \code{number}
#' @param gr The \code{\link{GRanges}} object
#' @param chunks The total number of desired chunks (rounded up to nearest int)
#' @param percent The percentage of overall CpGs in each chunk
#' @param size The number of CpGs to include in each chunk
#' @return \code{\link{GRangesList}} containing all the chunked \code{\link{GRanges}}
#' @import GenomicRanges
#' @examples
#' regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,100))
#' regions <- bin_granges(regions,bin_size=1)
#' split_granges(regions, chunks=10)
#' split_granges(regions, percent=10)
#' split_granges(regions, size=10)
#' @export
split_granges = function(gr,chunks = NA, percent = NA, size = NA) { #=NULL, percent = NULL

  if (!is(gr, "GRanges")) stop("Input must be a Granges object")
  
  if (sum(is.na(c(chunks,percent,size))) != 2) stop("Invalid input. Must contain 1 of either chunks, percent, or size")
  
  if (!is.na(percent)) {
    size <- floor(length(gr)*percent/(100))
  }
  
  if (!is.na(chunks)) {
    chunks = ceiling(chunks)
    size <- ceiling(length(gr)/chunks)
  }
  
  if (!is.na(size)) {
    splits <- ceiling(length(gr)/size)
    splits <- size*(0:(splits-1))+1
  }
  
  grl <- list()
  
  for (i in 1:length(splits)) {
    s <- splits[i]
    grl[[i]] <- gr[s:min(s+size-1,length(gr))]
  }
    
  return (GenomicRanges::GRangesList(grl))
}

#' Bins each region in a \code{\link{GRanges}} object into bins of specified \code{bin_size} 
#' @details Bins a single region in \code{\link{GRanges}} format into multiple regions with a specified \code{bin_size}. If \code{length(gr) %% bin_size != 0}, then the last GRange will have a length < \code{bin_size}. This is used instead of tile when you need consistently sized bins with the last bin being smaller
#' @param gr The \code{\link{GRanges}} object
#' @param bin_size The length of region in each bin
#' @return \code{\link{GRangesList}} containing all the binned \code{\link{GRanges}}
#' @import GenomicRanges
#' @examples
#' regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,10000))
#' bin_granges(regions,bin_size=1000) 
#' @export
bin_granges <- function(gr, bin_size = 100000) {#, enforce_size = FALSE) {

  if (!is(gr, "GRanges")) stop("Input must be a Granges object")
  
  ends <- len <- seqnames(gr)@lengths
  for (i in 1:length(ends)) ends[i] <- sum(as.vector(len[1:i]))
  starts <- head(c(1, ends + 1), -1)
  
  rngs <- lapply(1:length(starts), function (i) {
    
    #Get span of each chr     
    rng <- c(gr[starts[i]],gr[ends[i]])
    rng <- c(rng,gaps(rng,start=gr[starts[i]]@ranges@start))
    rng <- reduce(rng)
    
    #Create bins
    rng <- tile(rng, width = bin_size)
    unlist(rng)
  })
  
  rngs <- unlist(as(rngs, "GRangesList"))
  return(rngs)
  
  #TODO: implement enforce_size to make the minimum IRange for each chr be the length(bin_size)
}

#' Casts genomic regions into \code{\link{GRanges}} format
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
  unlockBinding("time.all",env=topenv()) #TODO: This looks hacky. Fix later
  unlockBinding("time.split",env=topenv())
  assign("time.all", proc.time()["elapsed"], envir=topenv())
  assign("time.split", proc.time()["elapsed"], envir=topenv())
  invisible(NULL)
}

#' Outputs the split/lap/iteration time
#' @details Gets the stored elapsed \\code{\link{proc.time}} from either the initial 
#' \code{\link{start_time}} or the previous \code{split_time}
#' @return Returns formatted elapsed time since \code{\link{start_time}} or last \code{\link{split_time}}
#' @export
split_time <- function() {
  time <- get("time.split", envir=topenv())
  if (!is.numeric(time)) {
    warning("start_time() not set. Starting from now.")
    start_time()
    return("[unknown time]")
  }
  time <- proc.time()["elapsed"]-time
  assign("time.split", proc.time()["elapsed"], envir=topenv())
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}

#' Stops an internal stopwatch and outputs overall time
#' @details Gets the stored elapsed \code{proc.time()} from initial \code{\link{start_time}} to calculate 
#' overall runtime
#' @return Returns formatted elapsed time since \code{\link{start_time}}
#' @export
stop_time <- function() {
  time <- get("time.all", envir=topenv())
  if (!is.numeric(time)) {
    warning("start_time() not set")
    return("[unknown time]")
  }
  time <- proc.time()["elapsed"]-get("time.all", envir=topenv())
  assign("time.split", NA, envir=topenv())
  assign("time.all", NA, envir=topenv())
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}

#' Subsets a given list of CpGs by another list of CpGs
#' @details Typically used to reduce the number of potential CpG sites to include only those present  in the input files so as to maximize performance and minimize resources. Can also be used for quality control to see if there is excessive number of CpG sites that are not present in the reference genome.
#' @param ref_cpgs A reference set of CpG sites (e.g. Hg19 or mm10) in bedgraph format
#' @param gen_cpgs A subset of CpG sites. Usually obtained from \code{\link{read_index}}.
#' @param verbose flag to output messages or not
#' @return Returns list of CpG sites in bedgraph format
#' @examples
#' ref_cpgs = data.frame(chr="chr1",start=(1:5*2-1), end=(1:5*2))
#' subset_ref_cpgs(ref_cpgs,ref_cpgs[1:3,])
#' @export
subset_ref_cpgs <- function(ref_cpgs, gen_cpgs, verbose = TRUE) {
  
  id <- NULL
  
  keys <- rbind(ref_cpgs[,c("chr","start")], gen_cpgs[,c("chr","start")])
  setDT(keys)[, id := .GRP, by = c("chr","start")]
  
  ref <- nrow(ref_cpgs)
  gen <- nrow(gen_cpgs)

  keys <- list(
    ref = keys[seq_len(ref),id],
    sub = keys[ref + seq_len(gen),id]
  )
  
  sub_cpgs <- ref_cpgs[keys$ref %in% keys$sub, , drop = FALSE]
  sub <- nrow(sub_cpgs)
  
  if (verbose) message("Dropped ",ref-sub,"/",ref," CpGs (",round((ref-sub)/ref*100,2),"%) from the reference set")
  if (verbose) message(gen-sub,"/",gen," subset CpGs (",round((gen-sub)/gen*100,2),"%) were not present in the reference set")
  
  return(sub_cpgs)
}

#' Gets bedgraph column indexes from common pipeline output formats
#' @details Typically used to reduce the number of potential CpG sites to include only those present  in the input files so as to maximize performance and minimize resources. Can also be used for quality control to see if there is excessive number of CpG sites that are not present in the reference genome.
#' @param protocol string; the protocol used for bedgraph output. Options are: "Bismark_cov", "MethylDackel", "MethylcTools", "BisSNP", "BSseeker2_CGmap"
#' @return List of column names and indexes
#' @export
get_source_idx = function(protocol = NULL) {
  
  protocol <- match.arg(
    arg = protocol,
    choices = c("Bismark_cov", "MethylDackel", "MethylcTools", "BisSNP", "BSseeker2_CGmap")
  )
  
  if (protocol == "MethylcTools") {
    return(list(col_idx = list(character = 1, numeric = 2, numeric = c(6, 7)),
                col_names = c("chr", "start", "M", "U"),
                fix_missing = c("cov := M+U", "beta := M/cov"), select = TRUE))
  } else if(protocol ==  "BisSNP"){
    return(list(col_idx = list(character = 1, numeric = 2, numeric = 4, numeric = 5),
                col_names = c("chr", "start", "beta", "cov"),
                fix_missing = c("M := as.integer(cov * beta)",
                                "U := cov - M"), select = TRUE))
  }else if(protocol ==  "BSseeker2"){
    return(list(col_idx = list(character = 1, numeric = 3, character = 4, numeric = 6, numeric = 8),
                col_names = c("chr", "start", "context", "beta", "cov"),
                fix_missing = c("context %in% 'CG'","M := as.integer(cov * beta)",
                                "U := cov - M"), select = TRUE))
  } else {
    # Bismark and methyldackel have same output format
    return(list(col_idx = list(character = 1, numeric = 2, numeric = 4, numeric = 5, numeric = 6),
                col_names = c("chr", "start", "beta", "M", "U"),
                fix_missing = c("cov := M+U"), select= TRUE))
  }
}

#' Subsets a given list of CpGs by another list of CpGs
#' @details Typically used to reduce the number of potential CpG sites to include only those present  in the input files so as to maximize performance and minimize resources. Can also be used for quality control to see if there is excessive number of CpG sites that are not present in the reference genome.
#' @param chr_idx integer; column of the chromosome
#' @param start_idx integer; column of the CpG start site
#' @param end_idx integer; column of the CpG end site
#' @param strand_idx integer; column of the strand
#' @param beta_idx integer; column of the beta value
#' @param M_idx integer; column of the # of methylated reads
#' @param U_idx integer; column of the # of unmethylated reads
#' @param cov_idx integer; column of the coverage
#' @param verbose flag to output messages or not
#' @return List of column names and indexes
#' @export
parse_source_idx = function(chr_idx = NULL, start_idx = NULL, end_idx = NULL, strand_idx = NULL,
                            beta_idx = NULL, M_idx = NULL, U_idx = NULL,
                            cov_idx = NULL, verbose = TRUE) {
  
  # mandatory chr and start field
  if (is.null(chr_idx) | is.null(start_idx)) {
    stop("missing chromosome/start indices\nUse pipeline argument if the files are from Bismark, MethyDeckal, or MethylcTools",
         call. = FALSE)
  }
  
  # See if any indices are duplicated
  if (length(which(duplicated(c(chr_idx, start_idx, end_idx, strand_idx, beta_idx, M_idx,
                                U_idx, cov_idx)))) > 0) {
    stop("Duplicated indices.", call. = FALSE)
  }
  
  # Check maximum betavalues (Can be 1 or 100)
  fix_missing = vector()
  if (is.null(strand_idx)) {
    fix_missing = "strand := '*'"
  }
  
  has_cov <- TRUE
  
  has <- function(x) {!is.null(x)}

  if (has(beta_idx)) {
    if (has(cov_idx)) {
      if (verbose) message("Estimating M and U from beta and cov")
      fix_missing = c(fix_missing,"M := as.integer(cov * beta)",
                      "U := cov - M")
    } else {
      if(has(M_idx) && has(U_idx)) { 
        if (verbose) message("Estimating cov from M and U")
        fix_missing = c(fix_missing, "cov := M+U") # Has: beta,M,U   No: cov
      } else if (has(M_idx)) {
        if (verbose) message("Estimating cov and U from M and beta")
        fix_missing = c(fix_missing, "cov := as.integer(M/beta)", # Has: beta,M,U   No: cov
                        "U := cov-M")
      } else if (has(U_idx)) {
        if (verbose) message("Estimating cov and M from U and beta")
        fix_missing = c(fix_missing, "cov := as.integer(U*(1-beta))", # Has: beta,U   No: M,cov
                        "M := cov-U")
      } else {
        if (verbose) message("Only beta info found")
        fix_missing = c(fix_missing, "cov := 2",
                        "M := beta*cov","U := (1-beta)*cov") # Has: beta   No: cov,M,U
        has_cov = FALSE
      }
    }
  } else { 
    if (has(cov_idx)) {
      if (has(M_idx)) {
        if (verbose) message("Estimating beta and U from M and cov. Beta will be [0,1]")
        fix_missing = c(fix_missing, "U := cov-M", "beta := M/cov")
      } else if (has(U_idx)) {
        if (verbose) message("Estimating beta and M from U and cov. Beta will be [0,1]")
        fix_missing = c(fix_missing, "M := cov-U", "beta := M/cov")
      } else {
        stop("Missing beta and cannot derive due to missing M and U.", call. = FALSE)
      }
    } else {
      if (is.null(M_idx) || is.null(U_idx)) {
        stop("Missing beta and cannot derive due to missing cov and one of M or U.", call. = FALSE)
      } else {
        if (verbose) message("Estimating beta and cov from M and U. Beta will be [0,1]")
        fix_missing = c(fix_missing, "cov := M+U", "beta := M/cov")
      }
    }
  }
   
  return(list(col_idx = c(chr = chr_idx, start = start_idx, end = end_idx, 
                            strand = strand_idx, beta = beta_idx, M = M_idx, 
                            U = U_idx,                                       
                            cov = cov_idx),
                fix_missing = fix_missing,
                has_cov = has_cov,
                select = FALSE))
}