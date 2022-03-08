#---- .validateArg -------------------------------------------------------------------------------------------
#' Validates inputs with multiple valid arguments. Allows partial matching.
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
#' If the argument for values is acceptable (e.g. "apple" or "ban"), it will return the matched string. If there is multiple matches (e.g., 'file' could match to c('file','filetype','rawfile')), it will still match if the strings exactly match.
#' 
#' This can be used with pipes, but will give erroneous names for the input variable
#' 
#' @param choices function definition or list of strings; This can be a function definition, and the argument with the identical name as the inputted 'arg' variable will be used, or can be a list of strings. If NULL, this will automatically assume that the parent function for which this function is called will be the function definition.
#' @param arg variable; the variable in which to check
#' @param ignore.case boolean; ignores case of the choices
#' @param partial.match boolean; whether to allow partial matching
#' @param multiple.match boolean; when to allow multiple input args. By default (FALSE), only the first input arg is used.
#' @return arg, if the value is in the function definition.
#' @export
.validateArg <- function(arg, choices = NULL, ignore.case = TRUE, partial.match = TRUE, multiple.match = FALSE) {

  #---- Input validation ---------------------------------------------------
  #.validateType(ignore.case,"boolean")
  #.validateType(partial.match,"boolean")
  #.validateType(multiple.match,"boolean")

  if (is.null(choices)) {
    choices <- deparse(sys.calls()[[sys.nframe()-1]])
    choices <- unlist(strsplit(choices, "[(]"))[[1]]
  } else if (is.function(choices)) {
    name = substitute(arg)
    choices <- eval(formals(choices)[[name]])
  }
  
  #---- Function code ------------------------------------------------------
  
  if (!multiple.match) {
    arg <- head(arg,1)
  } 

  if (any(sapply(choices, is.na)) && is.na(arg)) return (NA)

  matches <- NULL
  
  for (n in 1:length(arg)) {
  
    if (partial.match) {
      if (ignore.case) {
        m <- grepl(tolower(arg[n]), tolower(choices), fixed = TRUE)
      } else {
        m <- grepl(arg[n], choices, fixed = TRUE)
      }
    } else {
      if (ignore.case) {
        m <- choices %in% arg[n]
      } else {
        m <- tolower(choices) %in% tolower(arg[n])
      }
    }
    
    if (sum(m) != 1) stop(paste0("Invalid arg input for '",paste(name),"'. Found: '",arg[n],"'; Must match: '",
                                 paste0(choices, collapse="', '"),"'"), call. = FALSE)

    matches <- append(matches, choices[which(m)])
    
  }
    
  return(matches)
}

#---- .validateAssay -----------------------------------------------------------------------------------------
#' Validates an assay is in an [scMethrix()] experiment, using partial pattern matching.
#' 
#' This function checks whether an assay is present or absent within an experiment object, and will throw either an error or return a boolean based on results.
#' @param scm [scMethrix()]; the experiment object
#' @param assay string; the name of the assay
#' @param is.absent boolean; flag to check if the assay is absent (is.absent=TRUE), instead of present (is.absent=FALSE). Default = FALSE.
#' @param throws boolean; If validate fails, throw an error (throws=TRUE) or return FALSE (throws=FALSE). Default = TRUE.
#' @return string or boolean; if \code{check.absent == T}, the name of the matched assay or error if it doesn't exist. If \code{check_assay == F}, the boolean value for if the assay exists in the experiment
#' @export
.validateAssay <- function(scm = NULL, assay = NULL, is.absent = FALSE, throws = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  .validateType(assay,"string")
  .validateType(is.absent, "boolean")
  
  #---- Function code ------------------------------------------------------  
  if (!is.absent) {
    assay <- tryCatch(
      match.arg(arg = assay, choices = SummarizedExperiment::assayNames(scm)),
      error=function(cond) {
        if (throws) {
          stop(paste0("Invalid assay. No assay named '",assay,"' found in the experiment"), call. = FALSE)
        } else {
          return (invisible(FALSE))
        }
      }
    )
    return(invisible(assay))
  } else {
    return(!(assay %in% SummarizedExperiment::assayNames(scm)))
  }
}

#---- .validateType ------------------------------------------------------------------------------------------
#' Validates an input type
#' @param input the variable(s) to check. List structures are accepted and will be parsed. The input must be stored in a variable before passing to the function.
#' @param type [list()] of strings; the types to check for.
#' @param throws boolean; flag should the function throw an error (throws=TRUE) or return FALSE (throws=FALSE) if the validate fails. Default TRUE.
#' @param recursive_sub string; This is an internal variable used when recursively looping through lists. Do not modify the value.
#' @return invisible(TRUE) if validation passes. Error if validation fails (with throws=T) or invisible(FALSE) (with throws=T)
#' @export
#' @examples
#' testvar = "test"
#' print(.validateType(testvar,"string"))          # TRUE
#'
#' testvar = 5
#' tryCatch({
#'   .validateType(testvar,"string")               # ERROR
#' }, error = function(e) {
#'   print("Error!")
#' })
#' 
#' testvar = 5
#' print(.validateType(testvar,"string",throws=FALSE)) # FALSE
#'
#' testvar = list("test","test2")
#' print(.validateType(testvar,"string","list"))   # TRUE
.validateType <- function(input = NULL, type=c("Integer","Numeric","Character","String","Boolean","Logical","Vector",
                                               "List","File","Directory","GRanges","GenomicRanges","Function","Null",
                                               "NA","Dataframe","DF","S4","Distance","Chain","Soft"), throws=T, recursive_sub = NULL) {

  #---- Input validation ---------------------------------------------------
  if (length(type) == length(eval(formals(.validateType)[["type"]]))) {
    stop("No valid type specified.")
  }
  
  types <- sapply(type,function(type) .validateArg(type,.validateType))
  
  if (is.null(recursive_sub)) recursive_sub = gsub('"', "'", deparse(substitute(input)))
  
  #---- Function code ------------------------------------------------------
  valid <- F

    # Check for list structures
  if (length(input) > 1) {
    for (type in types) {
      if(type == "List"){
        valid <- is.list(input)
      } else if(type == "GRanges" | type == "GenomicRanges"){
        valid <- is(input, "GRanges")
      } else if (type == "Dataframe" | type == "DF") {
        valid <- is.data.frame(input)
      } else if (type == "Distance") {
        valid <- is(input,"dist")
      } else if (type == "Chain") {
        valid <- is(input,"Chain")
      }
    }
    if (valid) return(invisible(TRUE))
  }
  
  if (length(input) <= 1) {
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
      } else if(type == "List"){
        valid <- is.list(input)
      } else if (type == "File") {
        valid <- utils::file_test("-f", input)
      } else if (type == "Directory") {
        valid <- utils::file_test("-d", input)
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
      } else if (type == "Chain") {
          valid <- is(input,"Chain")
      } else if (type == "Soft") {
        valid <- (class(input) == "GSE")
      } else {
        stop("Invalid type with '",type,"'. This type is not supported for validation.", call. = FALSE)
      }
      
      if (valid) {
        break
      } else if (type == types[length(types)]) {
        if (throws) {
          stop("Invalid type input for '",recursive_sub,"'. Must be of type: '",
               paste0(types, collapse="', '"),"'", call. = FALSE)
        } else {
          return(invisible(FALSE)) 
        }
      }
    } 
  } else {valid = any(sapply(input, .validateType, type = types, throws = throws, recursive_sub = recursive_sub))}
  
  return(invisible(valid))
}

#---- .validateExp -------------------------------------------------------------------------------------------
#' Validates to see if object is a proper scMethrix object
#' @param scm [scMethrix()]; the experiment object to test
#' @param h5_dir string; the directory of an experiment object
#' @param throws boolean; whether to throw an error on a missing experiment. Will return FALSE on missing otherwise.
#' @return invisible(TRUE), if the object is valid. Error or FALSE if not.
#' @export
.validateExp <- function(scm = NULL, h5_dir = NULL, throws = TRUE) {
  
  if (!is.null(h5_dir)) {
    
    .validateType(h5_dir,"directory")

    tryCatch(
      expr = {
        scm <- load_scMethrix(h5_dir)
      },
      error = function(e){ 
        if (throws) {stop("Invalid scMethrix object found at ",h5_dir, call. = FALSE)
        } else {return(invisible(FALSE))}
      }
    )
  }
  
  if (!is.null(scm)) {
    if (!is(scm, "scMethrix")) {
      if (throws) {stop(paste0("Invalid scMethrix object supplied for '",substitute(scm),"'"), call. = FALSE)
      } else {return(invisible(FALSE))}
    }
  }
  
  return(invisible(TRUE))
}

#---- .validateValue -----------------------------------------------------------------------------------------
#' Validates numeric values based on some experession
#' @param value numeric; the value to test
#' @param ... string; the expressions to test
#' @return invisible(TRUE), if the object is valid. Error if not.
#' @export
.validateValue <- function(value,...) {
  
  if (!is.null(value) && !is.na(value)) {
    
    if (!.validateType(value,c("numeric","integer"),throws = FALSE))
      stop("Invalid value for '",substitute(value),"'. Must be of 'numeric' or 'integer' type.")
    
    conds = c()

    for (condition in list(...)) {
      if (!(eval(parse(text=paste0(value,condition))))) {
        conds <- c(conds,condition)
      }
    }
    
    if (length(conds != 0)) {
      conds <- paste(substitute(value),conds, collapse="', '")
      stop ("Invalid value: '",substitute(value)," = ",value,"'. Must fit condition: '",conds,"'", call. = FALSE)
    }
    
  }
  return(invisible(TRUE))
}

#---- .validateThreads ---------------------------------------------------------------------------------------
#' Validates the number of threads for the session. Windows can only support one thread
#' @param n_threads numeric; the number of threads
#' @return integer; 1 if windows, or some number of threads between 1 and parallel::detectCores
.validateThreads <- function(n_threads) {
  
  .validateType(n_threads,"integer")
  
  # if (grepl("Windows", Sys.getenv("OS"))) {
  #   if (n_threads > 1) warning("Invalid threads. Parallel processing is not enabled for non-POSIX system (e.g., Windows). Defaulting to 1.")
  #   return(1)
  # }

  return(max(min(parallel::detectCores(),n_threads),1))
}



#' Validates that a package is installed
#' @param package string; the name of the package
#' @return `ERROR` if not installed, `invisible(TRUE)` if installed.
.validatePackageInstall <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) 
    stop("Package '",package,"' is necessary for this function. Use: install.packages('",
         package,"')", call. = FALSE)
  
  return(invisible(TRUE))
}


#' Empty Value
#'
#' Rails-inspired helper that checks if vector values are "empty", i.e. if it's: \code{NULL}, zero-length, \code{NA}, \code{NaN}, \code{FALSE}, an empty string or \code{0}. Note that unlike its native R \code{is.<something>} sibling functions, \code{is.empty} is vectorised (hence the "values").
#' 
#' @references Take from Rapporter\\rapportools <https://github.com/Rapporter/rapportools> 
#' @param x an object to check its emptiness
#' @param trim trim whitespace? (\code{TRUE} by default)
#' @param ... additional arguments for \code{\link{sapply}}
#' @examples \dontrun{
#' is.empty(NULL)     # [1] TRUE
#' is.empty(c())      # [1] TRUE
#' is.empty(NA)       # [1] TRUE
#' is.empty(NaN)      # [1] TRUE
#' is.empty("")       # [1] TRUE
#' is.empty(0)        # [1] TRUE
#' is.empty(0.00)     # [1] TRUE
#' is.empty("    ")   # [1] TRUE
#' is.empty("foobar") # [1] FALSE
#' is.empty("    ", trim = FALSE)    # [1] FALSE
#' # is.empty is vectorised!
#' all(is.empty(rep("", 10)))        # [1] TRUE
#' all(is.empty(matrix(NA, 10, 10))) # [1] TRUE
#' }
#' @export
is.empty <- function(x, trim = TRUE, ...) {
  if (length(x) <= 1) {
    if (is.null(x))
      return (TRUE)
    if (length(x) == 0)
      return (TRUE)
    if (is.na(x) || is.nan(x))
      return (TRUE)
    if (is.character(x) && nchar(ifelse(trim, trim.space(x), x)) == 0)
      return (TRUE)
    if (is.logical(x) && !isTRUE(x))
      return (TRUE)
    if (is.numeric(x) && x == 0)
      return (TRUE)
    return (FALSE)
  } else
    sapply(x, is.empty, trim = trim, ...)
}

#' Trim Spaces
#'
#' Removes leading and/or trailing space(s) from a character vector. By default, it removes both leading and trailing spaces.
#' 
#' @references Take from Rapporter\\rapportools <https://github.com/Rapporter/rapportools> 
#' @param x a character vector which values need whitespace trimming
#' @param what which part of the string should be trimmed. Defaults to \code{both} which removes trailing and leading spaces. If \code{none}, no trimming will be performed.
#' @param space.regex a character value containing a regex that defines a space character
#' @param ... additional arguments for \code{\link{gsub}} function
#' @return a character vector with (hopefully) trimmed spaces
trim.space <- function(x, what = c('both', 'leading', 'trailing', 'none'), space.regex = '[:space:]', ...){
  if (missing(x))
    stop('nothing to trim spaces to =(')
  re <- switch(match.arg(what),
               both     = sprintf('^[%s]+|[%s]+$', space.regex, space.regex),
               leading  = sprintf('^[%s]+', space.regex),
               trailing = sprintf('[%s]+$', space.regex),
               none     = {
                 return (x)
               })
  vgsub(re, '', x, ...)
}

#' Vectorised String Replacement
#'
#' A simple wrapper for \code{\link{gsub}} that replaces all patterns from \code{pattern} argument with ones in \code{replacement} over vector provided in argument \code{x}.
#' @references See original thread for more details \url{http://stackoverflow.com/a/6954308/457898}. Special thanks to user Jean-Robert for this one!
#' @param pattern see eponymous argument for \code{\link{gsub}} function
#' @param replacement see eponymous argument for \code{\link{gsub}} function
#' @param x see eponymous argument for \code{\link{gsub}} function
#' @param ... additional arguments for \code{\link{gsub}} function
#' @return a character vector with string replacements
vgsub <- function(pattern, replacement, x, ...){
  for(i in 1:length(pattern))
    x <- gsub(pattern[i], replacement[i], x, ...)
  x
}
