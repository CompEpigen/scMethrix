
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

#' Validates numeric values based on some experession
#' @param value numeric; the value to test
#' @param ... string; the expressions to test
#' @return invisible(TRUE), if the object is valid. Error if not.
.validateValue <- function(value,...) {
  
  if (!is.null(value) && !is.na(value)) {
    
    if (!.validateType(value,c("numeric","integer"),throws=F))
      stop("Invalid value for '",substitute(value),"'. Must be of 'numeric' or 'integer' type.")
    
    for (condition in list(...)) {
      if (!(eval(parse(text=paste0(value,condition))))) {
        stop ("Invalid value: '",substitute(value)," = ",value,"'. Must fit condition: ",substitute(value)," ",condition, call. = FALSE)
      }
    }
  }
  return(invisible(TRUE))
}

#' Validates the number of threads for the session. Windows can only support one thread
#' @param n_threads numeric; the number of threads
#' @return integer; 1 if windows, or some number of threads between 1 and parallel::detectCores
.validateThreads <- function(n_threads) {
  
  .validateType(n_threads,"integer")
  
  # if (grepl("Windows", Sys.getenv("OS"))) {
  #   if (n_threads > 1) warning("Invalid threads. Parallel processing is not enabled for non-POSIX system (e.g., Windows). ")
  #   return(0)
  # } 
  
  return(max(min(parallel::detectCores(),n_threads),1))
}