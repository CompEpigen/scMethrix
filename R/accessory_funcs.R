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

#' Bins each region in a \code{\link{GRanges}} object into bins of specified \code{bin_size} 
#' @details Bins a single region in \code{\link{GRanges}} format into multiple regions with a specified \code{bin_size}. If \code{length(gr) %% bin_size != 0}, then the last GRange will have a length < \code{bin_size}. This is used instead of tile when you need consistently sized bins with the last bin being smaller
#' @param gr GRanges; The \code{\link{GRanges}} object
#' @param bin_size integer; x > 0; The length of region in each bin
#' @return \code{\link{GRangesList}} containing all the binned \code{\link{GRanges}}
#' @import GenomicRanges
#' @examples
#' regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,10000))
#' bin_granges(regions,bin_size=1000) 
#' @export
bin_granges <- function(gr, bin_size = 100000) {#, enforce_size = FALSE) {

  #- Input Validation --------------------------------------------------------------------------
  .validateType(gr,"Granges")
  .validateType(bin_size, "integer")
  .validateValue(bin_size,">0")
  
  #- Function code -----------------------------------------------------------------------------
  gr <- GenomicRanges::slidingWindows(gr,width=bin_size,step=bin_size)
  return(unlist(as(gr, "GRangesList")))
}

#' Casts genomic regions into \code{\link{GRanges}} format
#' @details Casts the input as a \code{\link{GRanges}} object. Input can be \code{\link{GRanges}} or a 
#' \code{\link{data.frame}}-compatible class that can be cast through \code{as.data.frame()}. Input BED format
#'  must be \code{chr-start-end} for \code{\link{data.frame}} objects.
#' @param regions GRanges or data.frame; The input regions to cast to \code{\link{GRanges}}
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
  assign("time.all", proc.time()["elapsed"], envir=timer.env)
  assign("time.split", proc.time()["elapsed"], envir=timer.env)
  invisible(NULL)
}

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

#' Subsets a given list of CpGs by another list of CpGs
#' @details Typically used to reduce the number of potential CpG sites to include only those present  in the input files so as to maximize performance and minimize resources. Can also be used for quality control to see if there is excessive number of CpG sites that are not present in the reference genome.
#' @param ref_cpgs data.table; A reference set of CpG sites (e.g. Hg19 or mm10) in bedgraph format
#' @param gen_cpgs data.table; A subset of CpG sites. Usually obtained from \code{\link{read_index}}.
#' @param verbose boolean; flag to output messages or not
#' @return Returns list of CpG sites in bedgraph format
#' @examples
#' ref_cpgs = data.frame(chr="chr1",start=(1:5*2-1), end=(1:5*2))
#' subset_ref_cpgs(ref_cpgs,ref_cpgs[1:3,])
#' @export
subset_ref_cpgs <- function(ref_cpgs, gen_cpgs, verbose = TRUE) {
  
  #- data.table initialization -----------------------------------------------------------------
  id <- NULL
  
  #- Function code -----------------------------------------------------------------------------
  keys <- rbind(ref_cpgs[,c("chr","start")], gen_cpgs[,c("chr","start")])
  data.table::setDT(keys)[, id := .GRP, by = c("chr","start")]
  
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

#' Validates arguments. Allows partial matching.
#' @details Check the parent function input arguments to see whether the inputted value is part of the set. Will return a formatted error message with the incorrect variable name and all the acceptable inputs.
#' 
#' To use, it can be called as such:
#' 
#' genericFunc <- function(values = c("apple","orange","banana")) {
#' 
#'    values <- .validateArg(values)
#' 
#' }
#' 
#' If the argument for values is acceptable (e.g. "apple" or "ban"), it will return the matched string.
#' 
#' This can be used with pipes, but will give erroneous names for the input variable
#' 
#' @param parent closure; the parent function in which to check the input
#' @param arg variable; the variable in which to check
#' @param ignore.case boolean; ignores case of the choices
#' @return arg, if the value is in the function definition.
.validateArg <- function(arg, parent = NULL, ignore.case = T) {

  #.validateType(ignore.case,"boolean")
  
  #- Function code -----------------------------------------------------------------------------
  if (is.null(parent)) {
    parent <- deparse(sys.calls()[[sys.nframe()-1]])
    parent <- unlist(strsplit(parent, "[(]"))[[1]]
  }
  
  name = substitute(arg)
  choices <- eval(formals(parent)[[name]])
  
  arg <- tryCatch(
    {
      if (ignore.case) {
        m <- match.arg(arg = tolower(arg), choices = tolower(choices))
        i <- which(tolower(choices) %in% m)
        choices[i]
        
      } else {
        match.arg(arg = arg, choices = choices)
      }
    },
    error=function(cond) {
      stop(paste0("Invalid arg input for '",paste(name),"'. Found: '",arg,"'; Must match: '",
                  paste0(eval(formals(parent)[[name]]), collapse="', '"),"'"), call. = FALSE)
    }
  )
  
  return(arg)
}

#' Validates an assay is in the object. Allows partial pattern.
#' @details Check the assays in an scMethrix object and partial matches 
#' @param scm scMethrix; the experiment object
#' @param assay string; the name of the assay
#' @param check.absent boolean; Checks if the assay is present
#' @return string or boolean; if \code{check.absent == T}, the name of the matched assay or error if it doesn't exist. If \code{check_assay == F}, the boolean value for if the assay exists in the experiment
.validateAssay <- function(scm = NULL,assay = NULL, check.absent = F) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  .validateType(assay,"string")
  .validateType(check.absent, "boolean")

  #- Function code -----------------------------------------------------------------------------  
  if (!check.absent) {
    assay <- tryCatch(
      match.arg(arg = assay, choices = SummarizedExperiment::assayNames(scm)),
      error=function(cond) 
          stop(paste0("Invalid assay. No assay named '",assay,"' found in the experiment '",substitute(scm),"'"), call. = FALSE)
    )
    return(invisible(assay))
  } else {
    return(!(assay %in% SummarizedExperiment::assayNames(scm)))
  }
}

# .validateType <- function(input = NULL, type=c("Integer","Numeric","Character","String","Boolean","Logical","Vector",
#                                                "List","File","Directory","GRanges","GenomicRanges","Function","Null",
#                                                "Dataframe","DF","S4")) {
# 
#   if (b) browser()
#   
#   if (length(type) == length(eval(formals(.validateType)[["type"]]))) {
#     stop("No valid type specified.")
#   }
#   
#   types <- sapply(type,function(type) .validateArg(type,.validateType))
#   inputs <- input#list(input,NULL) # necessary to avoid iterating through iterable objects (e.g. GRanges)
# 
#   for (type in types) {
#     
#     #For container formats that are iterable
#     if (c("GRanges","GenomicRanges","S4")) {
#       
#       valid <- F
#       
#       for (type in types) {
#         
#         if(type == "Vector"){
#           valid <- is.vector(input)
#         } else if(type == "List"){
#           valid <- is.list(input)
#         } else if(type == "GRanges" | type == "GenomicRanges"){
#           valid <- is(input, "GRanges")
#         } else if (type == "Dataframe" | type == "DF") {
#           valid <- is.data.frame(input) 
#         } else if (type == "S4") {
#           valid <- isS4(input) 
#         } else {
#           stop("Invalid type with '",type,"'. This type is not supported. Also, this should never be reached.")
#         }
#         
#       }
#       
#       if (!valid) {
#         stop("Invalid type input for '",substitute(input),"'. Must be of type: '",
#              paste0(type, collapse="', '"),"'", call. = FALSE)
#       }
#     }
#       
#       
#     } else {
#       for (input in inputs) {#inputs[1:(length(inputs)-1)]) {
#         
#         valid <- F
#         
#         for (type in types) {
#           
#           if (type == "Integer") {
#             if(is.numeric(input)) valid <- (input == round(input))
#           } else if (type == "Numeric") {
#             valid = is.numeric(input) 
#           } else if (type == "Character") {
#             valid = is.character(input) && nchar(input)==1
#           } else if (type == "String") {
#             valid = is.character(input) 
#           } else if(type == "Boolean" | type == "Logical"){
#             valid <- is.logical(input)
#           } else if (type == "File") {
#             valid <- all(file.exists(input))
#           } else if (type == "Directory") {
#             valid <- all(dir.exists(input))
#           } else if (type == "Function") {
#             valid <- is.function(i)
#           } else if (type == "Null") {
#             valid <- is.null(input) 
#           } else {
#             stop("Invalid type with '",type,"'. This type is not supported. Also, this should never be reached.")
#           }
#           
#         }
#         
#         if (!valid) {
#           stop("Invalid type input for '",substitute(input),"'. Must be of type: '",
#                paste0(type, collapse="', '"),"'", call. = FALSE)
#         }
#       }
#       
#     }
#     
#     
#     
#     
#     
#   }
#   
#   
#   
#   
#   
#   
#   return(invisible(TRUE))
# }


.validateType <- function(input = NULL, type=c("Integer","Numeric","Character","String","Boolean","Logical","Vector",
                                               "List","File","Directory","GRanges","GenomicRanges","Function","Null",
                                               "NA","Dataframe","DF","S4","Distance"), throws=T) {
    
  #- Input Validation --------------------------------------------------------------------------
  if (length(type) == length(eval(formals(.validateType)[["type"]]))) {
    stop("No valid type specified.")
  }

  types <- sapply(type,function(type) .validateArg(type,.validateType))
  # inputs <- input#list(input,NULL) # necessary to avoid iterating through iterable objects (e.g. GRanges)s
  
  #- Function code -----------------------------------------------------------------------------
  valid <- F

  for (type in types) {
    if (type == "Null") {
      valid <- is.null(input)
    }  else if (type == "NA") {
      valid <- is.na(input)
    } else if (type == "Integer") {
      if(is.numeric(input)) valid <- (input == round(input))
    } else if (type == "Numeric") {
      valid = is.numeric(input)
    } else if (type == "Character") {
      valid = is.character(input) && nchar(input)==1
    } else if (type == "String") {
      valid = is.character(input)
    } else if(type == "Boolean" | type == "Logical"){
      valid <- is.logical(input)
    } else if(type == "Vector"){
      valid <- is.vector(input)
    } else if(type == "List"){
      valid <- is.list(input)
    } else if (type == "File") {
      valid <- all(file.exists(input))
    } else if (type == "Directory") {
      valid <- all(dir.exists(input))
    } else if(type == "GRanges" | type == "GenomicRanges"){
      valid <- is(input, "GRanges")
    } else if (type == "Function") {
      valid <- is.function(input)
    } else if (type == "Dataframe" | type == "DF") {
      valid <- is.data.frame(input)
    } else if (type == "S4") {
      valid <- isS4(input)
    } else if (type == "Distance") {
      valid <- is(input,"dist")
    } else {
      stop("Invalid type with '",type,"'. This type is not supported for validation.")
    }

    if (valid) {
      break
    } else if (type == types[length(types)]) {
      if (throws) {
      stop("Invalid type input for '",substitute(input),"'. Must be of type: '",
           paste0(types, collapse="', '"),"'", call. = FALSE)
      } else {
        return(invisible(FALSE)) 
      }
    }
  }

  return(invisible(TRUE))
}


#' Validates to see if object is a proper scMethrix object
#' @param scm scMethrix; the experiment object to test
#' @return invisible(TRUE), if the object is valid. Error if not.
.validateExp <- function(scm) {
  if (!is(scm, "scMethrix")) stop(paste0("Invalid scMethrix object supplied for '",substitute(scm),"'"), call. = FALSE)
  return(invisible(TRUE))
}

.validateValue <- function(value,...) {

    if (!is.null(value) && !is.na(value)) {
      
      if (!.validateType(value,c("numeric","integer"),throws=F))
        stop("Invalid value for '",substitute(value),"'. Must be of 'numeric' or 'integer' type.")
      
      for (condition in list(...)) {
        if (!(eval(parse(text=paste0(value,condition))))) {
          stop ("Invalid value for '",substitute(value),"'. Must fit condition: ",value,condition)
        }
      }
    }
  return(invisible(TRUE))
}