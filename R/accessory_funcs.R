#--- is_h5 --------------------------------------------------------------------------------------------------
#' Checks if \code{\link{scMethrix}} object is an HDF5 object
#' @param scm scMethrix; The \code{\link{scMethrix}} object
#' @return boolean Whether the object is HDF5
#' @examples
#' data('scMethrix_data')
#' is_h5(scMethrix_data)
#' @export
is_h5 = function(scm) {
  .validateExp(scm)
  return(scm@metadata$is_h5)
}

#--- has_cov ------------------------------------------------------------------------------------------------
#' Checks if \code{\link{scMethrix}} object has a coverage matrix
#' @param scm scMethrix; The \code{\link{scMethrix}} object
#' @return boolean Whether the object has a coverage matrix
#' @examples
#' data('scMethrix_data')
#' has_cov(scMethrix_data)
#' @export
has_cov = function(scm) {
  .validateExp(scm)
  return("counts" %in% SummarizedExperiment::assayNames(scm))
}

#--- get_sample_name ----------------------------------------------------------------------------------------
#' Returns sample name derived from the input file name
#' @param s string; A file.path
#' @return string containing the sample name
#' @import tools
#' @examples
#' #get_sample_name("C:/dir/dir/filename.ext")
#' @export
get_sample_name = function(s) {
  if (!is.character(s)) stop("Must be a string file path")
  return(tools::file_path_sans_ext(basename(s)))
}

#--- binarize -----------------------------------------------------------------------------------------------
#' Binarize an input value based on a \code{threshold}
#' @details Assigns a value of 0 or 1 based on being < or > the \code{thresdhold}, respectively.
#'  If \code{x} = \code{threshold}, \code{x} = 0. NA values are assigned as \code{rep_na}
#' @param x numeric; A vector to binarize
#' @param threshold numeric; The threshold for binarizing. Will default to the center number between max and min.
#' @param rep_na numeric; The value to replace missing values with. Default NA. 
#' @param verbose boolean; flag for whether to display threshold or not
#' @return vector; contains {0,1}, if below or above the threshold, or 'rep.na' if NA
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

#--- fill ---------------------------------------------------------------------------------------------------
#' Fills a vector with a specified \code{fill} value
#' @param x vector; A vector in which to fill the NA values
#' @param fill basic data type; Any value from one of R's basic data types (character, numeric, integer, logical, complex)
#' @return vector; Same values as input vector, but NA values are replaced with \code{fill} if above of below the threshold, or 'rep.na' if NA
#' @examples
#' vals <- c(0,0.25,0.5,0.75,1,NA)
#' fill(vals, fill=2)
#' @export
fill = function(x, fill = NA) {
  x[is.na(x)] <- fill
  return(x)
}

#--- colbind ------------------------------------------------------------------------------------------------
#' A faster version of cbind when trying to combine lists of data.tables
#' @param ... A list of data.tables with identical # of rows
#' @return data.table; the cbinded output 
#' @export
colbind = function(...) {
  setDT(
    unlist(..., recursive = FALSE),
    check.names = TRUE
  )[]
}

#--- split_vector -------------------------------------------------------------------------------------------
#' Splits a vector into list of vectors by \code{chunk} or \code{size}
#' @details Splits a vector into consistantly sized sub-lists. The sub-list size will always be decreasing based on list order.
#' @param vec vector; The vector to split
#' @param chunks integer; x > 0; The number of desired sub-lists
#' @param percent numeric; 100 > x > 0; The maximum percentage of elements each sub-list should hold
#' @param size integer; x > 0; The maximum size of each sub-list
#' @return A list of sub-vectors
#' @examples
#' # Split vector into 4 sub vectors
#' split_vector(c(1,2,3,4,5,6,7,8),chunks=4)
#' 
#' # Split vector into sub-vectors with a size of 2
#' split_vector(c(1,2,3,4,5,6,7,8),size=2)
#' 
#' # Split vector into sub-vectors each with 25% of the total elements
#' split_vector(c(1,2,3,4,5,6,7,8),percent=25)
#' @export
split_vector = function(vec, chunks = NA, percent = NA, size = NA) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateType(chunks,c("integer","na"))
  .validateType(size,c("integer","na"))
  .validateType(percent,c("numeric","na"))
  
  .validateValue(chunks,">0")
  .validateValue(size,">0")
  .validateValue(percent,">0","<100")
  
  if (sum(is.na(c(chunks,percent,size))) != 2) stop("Invalid input. Must contain 1 of either chunks, percent, or size")
  
  #- Function code -----------------------------------------------------------------------------  
  if (!is.na(percent)) chunks = 100/percent
  if (!is.na(size)) chunks = length(vec)/ceiling(size)
  chunks = max(1,chunks)
  return(unname(split(vec, sort(rep_len(1:ceiling(chunks), length(vec))))))
}

#--- start_time ---------------------------------------------------------------------------------------------
#' Starts an internal stopwatch
#' @details Save the current time to later use for split/lap and overall times
#' @return NULL
#' @export
start_time <- function() {
  assign("time.all", proc.time()["elapsed"], envir=timer.env)
  assign("time.split", proc.time()["elapsed"], envir=timer.env)
  invisible(NULL)
}

#--- split_time ---------------------------------------------------------------------------------------------
#' Outputs the split/lap/iteration time
#' @details Gets the stored elapsed \\code{\link{proc.time}} from either the initial
#' \code{\link{start_time}} or the previous \code{split_time}
#' @return Returns formatted elapsed time since \code{\link{start_time}} or last \code{\link{split_time}}
#' @export
split_time <- function() {
  time <- get("time.split", envir=timer.env)
  if (!is.numeric(time)) {
    warning("start_time() not set. Starting from now.")
    start_time()
    return("[unknown time]")
  }
  time <- proc.time()["elapsed"]-time
  assign("time.split", proc.time()["elapsed"], envir=timer.env)
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}

#--- stop_time ----------------------------------------------------------------------------------------------
#' Stops an internal stopwatch and outputs overall time
#' @details Gets the stored elapsed \code{proc.time()} from initial \code{\link{start_time}} to calculate
#' overall runtime
#' @return Returns formatted elapsed time since \code{\link{start_time}}
#' @export
stop_time <- function() {
  time <- get("time.all", envir=timer.env)
  if (!is.numeric(time)) {
    #warning("start_time() not set")
    return("[unknown time]")
  }
  time <- proc.time()["elapsed"]-get("time.all", envir=timer.env)
  assign("time.split", NA, envir=timer.env)
  assign("time.all", NA, envir=timer.env)
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}

#--- get_source_idx -----------------------------------------------------------------------------------------
#' Gets bedgraph column indexes from common pipeline output formats
#' @details Typically used to reduce the number of potential CpG sites to include only those present  in the input files so as to maximize performance and minimize resources. Can also be used for quality control to see if there is excessive number of CpG sites that are not present in the reference genome.
#' @param protocol string; the protocol used for bedgraph output. Options are: "Bismark_cov", "MethylDackel", "MethylcTools", "BisSNP", "BSseeker2_CGmap"
#' @return List of column names and indexes
#' @examples
#' get_source_idx("MethylDackel")
#' @export
get_source_idx = function(protocol = c("Bismark_cov", "MethylDackel", "MethylcTools", "BisSNP", "BSseeker2_CGmap")) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateArg(protocol, get_source_idx)
  
  #- Function code -----------------------------------------------------------------------------
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

#--- parse_source_idx ---------------------------------------------------------------------------------------
#' Generates the column structure for importing bedgraph files
#' @details 
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
  
  #- Input Validation --------------------------------------------------------------------------
  .validateType(chr_idx,   c("integer","null"))
  .validateType(start_idx, c("integer","null"))
  .validateType(end_idx,   c("integer","null"))
  .validateType(strand_idx,c("integer","null"))
  .validateType(beta_idx,  c("integer","null"))
  .validateType(M_idx,     c("integer","null"))
  .validateType(U_idx,     c("integer","null"))
  .validateType(cov_idx,   c("integer","null"))
  .validateType(verbose,   "boolean")
  
  .validateValue(chr_idx,   ">0")
  .validateValue(start_idx, ">0")
  .validateValue(end_idx,   ">0")
  .validateValue(strand_idx,">0")
  .validateValue(beta_idx,  ">0")
  .validateValue(M_idx,     ">0")
  .validateValue(U_idx,     ">0")
  .validateValue(cov_idx,   ">0")
  
  if (is.null(chr_idx) | is.null(start_idx)) {
    stop("Missing chromosome/start indices\nUse pipeline argument if the files are from Bismark, MethyDeckal, or MethylcTools",
         call. = FALSE)
  }
  
  if (length(which(duplicated(c(chr_idx, start_idx, end_idx, strand_idx, beta_idx, M_idx,
                                U_idx, cov_idx)))) > 0) {
    stop("Duplicated indices specified for the source bedgraph file.", call. = FALSE)
  }
  
  #- Function code -----------------------------------------------------------------------------
  fix_missing = vector()
  if (is.null(strand_idx)) {
    fix_missing = "strand := '*'"
  }
  
  has_cov <- TRUE
  
  has <- function(x) {!is.null(x)}
  
  if (has(beta_idx)) {
    if (has(cov_idx)) {
      if (verbose) message("   Estimating M and U from beta and cov")
      fix_missing = c(fix_missing,"M := as.integer(cov * beta)",
                      "U := cov - M")
    } else {
      if(has(M_idx) && has(U_idx)) { 
        if (verbose) message("   Estimating cov from M and U")
        fix_missing = c(fix_missing, "cov := M+U") # Has: beta,M,U   No: cov
      } else if (has(M_idx)) {
        if (verbose) message("   Estimating cov and U from M and beta")
        fix_missing = c(fix_missing, "cov := as.integer(M/beta)", # Has: beta,M,U   No: cov
                        "U := cov-M")
      } else if (has(U_idx)) {
        if (verbose) message("   Estimating cov and M from U and beta")
        fix_missing = c(fix_missing, "cov := as.integer(U*(1-beta))", # Has: beta,U   No: M,cov
                        "M := cov-U")
      } else {
        if (verbose) message("   Only beta info found")
        fix_missing = c(fix_missing, "cov := 2",
                        "M := beta*cov","U := (1-beta)*cov") # Has: beta   No: cov,M,U
        has_cov = FALSE
      }
    }
  } else { 
    if (has(cov_idx)) {
      if (has(M_idx)) {
        if (verbose) message("   Estimating beta and U from M and cov. Beta will be [0,1]")
        fix_missing = c(fix_missing, "U := cov-M", "beta := M/cov")
      } else if (has(U_idx)) {
        if (verbose) message("   Estimating beta and M from U and cov. Beta will be [0,1]")
        fix_missing = c(fix_missing, "M := cov-U", "beta := M/cov")
      } else {
        stop("Missing beta and cannot derive due to missing M and U.", call. = FALSE)
      }
    } else {
      if (is.null(M_idx) || is.null(U_idx)) {
        stop("Missing beta and cannot derive due to missing cov and one of M or U.", call. = FALSE)
      } else {
        if (verbose) message("   Estimating beta and cov from M and U. Beta will be [0,1]")
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
