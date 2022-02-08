#' Accessors methods and slot descriptions for [scMethrix()]
#' @name scMethrix-class
#' @description S4 class scMethrix
#' @docType class
#' @slot assays [list()]; assays containing methylation or coverage information
#' @slot colData data.frame; metadata corresponding to samples
#' @slot metadata [list()]; metadata pertaining to the experiment
#' @slot rowRanges [GenomicRanges::GRanges()]; the genomic coordinates of CpG sites and associated metadata
#' @exportClass scMethrix
#' @param object,x An [scMethrix()] object
#' @aliases 
#' coerce,GRset,scMethrix-method,scMethrix-method coerce
#' coerce,SingleCellExperiment,scMethrix-method,scMethrix-method coerce
#' @section Coercion:
#' An scMethrix object can be coerced from a [minfi::GenomicRatioSet()] or [SingleCellExperiment::SingleCellExperiment()] using the [methods::as()] function. However, due to limitations of [minfi::GenomicRatioSet()], coverage information cannot be recovered from a [minfi::GenomicRatioSet()]. As well, the output [scMethrix()] object will be in HDF5 format by default.
#' \preformatted{(x, "scMethrix")}
#' @importFrom stats median quantile sd density
#' @importFrom utils data head write.table menu browseURL
#' @importFrom methods is as new
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @seealso scMethrix
setClass(Class = "scMethrix", contains = "SingleCellExperiment")

#' [scMethrix()] - fast and easy summarization of methylation data into an integrative experiment object.
#' 
#' The object combines multiple data containers representing common data from methylation experiments (e.g., samples, features, assays). It inherits from the SingleCellExperiment class and is used in the same manner, but with additional consistency checks and methylation-specific functionality. 
#' @inheritParams generic_scMethrix_function
#' @param assays [list()] of matrices; The assays to include in the experiment
#' @param colData data.frame; The metadata corresponding to each sample
#' @param rowRanges [GenomicRanges::GRanges()]; The genomic loci corresponding to the assays
#' @param is_h5 boolean; Should assays be saved as HDF5 format?
#' @param genome string; The name of the genome build used. This will overwrite an element of the same name in \code{metadata}, if present. Default = NULL.
#' @param metadata named [list()] of strings; list of relevant experiment data. Elements with the name of 'is_h5' and 'genome' will be replaced by the arguments above. Default = NULL.
#' @export scMethrix
#' @seealso scMethrix-class for object structure and accessors
scMethrix <- function(assays = list(), colData = S4Vectors::DataFrame(), rowRanges = GenomicRanges::GRanges(), 
                      is_h5 = FALSE, h5_dir = NULL, genome = NULL, metadata = NULL, replace = FALSE, verbose=TRUE) {

  # Ensure consistent metadata
  if (is.null(genome) && "genome" %in% names(metadata)) genome <- metadata[["genome"]]
  metadata <- metadata[!names(metadata) == "is_h5"]
  metadata <- c(list(is_h5 = is_h5, genome = genome), metadata)
  
  # Save the experiment into memory or HDF5 format
  if (is_h5) {
    
    scm <- new("scMethrix",SingleCellExperiment(assays = lapply(assays,function(x) as(x, "HDF5Array")), 
                                                      colData = colData, rowRanges = rowRanges,
                                                      metadata = metadata))
    
    if (!is.null(h5_dir)) {
      tryCatch(save_scMethrix(scm = scm, dest = h5_dir, replace = replace, verbose = verbose),
               error = function(e) message(e,"\nThe dataset is not saved. Please save manually using the HDF5Array::saveSummarizedExperiment command."))
    }
    
  } else {
    scm <- new("scMethrix",SingleCellExperiment(assays = lapply(assays,function(x) as(x, "matrix")),
                                                      colData = colData, rowRanges = rowRanges,
                                                      metadata = metadata))
  }
  
  return(scm)
}

#---- Generic methods -------------------------------------------------------------------------------------------------
#' @describeIn scMethrix-class Show method for an [scMethrix()] object
setMethod(f = "show", signature = "scMethrix", definition = function(object) {
  
  h5 <- is_h5(object)
  if (h5) h5 <- path(score(object))
  
  cat(paste0("An object of class ", class(object), "\n"))
  cat(paste0("   CpGs: ", format(nrow(object), big.mark = ","),"\n"))
  cat(paste0("   Samples: ", ncol(object),"\n"))
  cat(paste0("   Assays: ", (paste(SummarizedExperiment::assayNames(object),collapse=", ")),"\n"))
  cat(paste0("   Reduced dims: ", (paste(SingleCellExperiment::reducedDimNames(object),collapse=", ")),"\n"))
  cat(paste0("   HDF5: ", h5, "\n"))
  cat(paste0("   Ref.Genome: ", S4Vectors::metadata(object)$genome, "\n"))
  cat(paste0("   Memory size: ", format(utils::object.size(object), units = "auto"), "\n"))
})

#---- showMore ---------------------------------------------------------------------------------------------------------
setGeneric("showMore", function(object) standardGeneric("showMore"))

#' @describeIn scMethrix-class Same as the show method, but adds additional information
setMethod(f = "showMore", signature = "scMethrix", definition = function(object)   {
  
  feature.names <- row.names(rowData(object))
  if (!is.null(feature.names)) feature.names <- paste0(substr(feature.names[1:5],1,8),"...",collapse=" l ")
  sample.names <- row.names(colData(object))
  if (!is.null(sample.names)) {
    sample.names <- substr(row.names(colData(object))[1:5],1,8)
    sample.names <- gsub('[^a-zA-Z]*$','',sample.names)
    sample.names <- paste0(sample.names,"...",collapse=" l ")
    sample.names <- paste0(" (",sample.names,")")
  }
  
  h5 <- is_h5(object)
  if (h5) h5 = path(score(object))
  
  cat(paste0("An object of class ", class(object), "\n"))
  cat(paste0("   CpGs: ", format(nrow(object), big.mark = ","),feature.names,"\n"))
  cat(paste0("   Samples: ", ncol(object),sample.names,"\n"))
  cat(paste0("   Assays: ", (paste(SummarizedExperiment::assayNames(object),collapse=", ")),"\n"))
  cat(paste0("   Reduced dims: ", (paste(SingleCellExperiment::reducedDimNames(object),collapse=", ")),"\n"))
  cat(paste0("   HDF5: ", h5, "\n"))
  cat(paste0("   Ref.Genome: ", S4Vectors::metadata(object)$genome, "\n"))
  cat(paste0("   Memory size: ", format(utils::object.size(object), units = "auto"), "\n"))
})

#---- score ------------------------------------------------------------------------------------------------------------
#' @describeIn scMethrix-class Gets the assay named "score"
setMethod(f = "score", signature = "scMethrix", definition = function(x) {get_matrix(x,"score")})

#---- counts -----------------------------------------------------------------------------------------------------------
#' @describeIn scMethrix-class Gets the assay named "counts"
setMethod(f = "counts", signature = "scMethrix", definition = function(object) {get_matrix(object,"counts")})

#---- has_cov ----------------------------------------------------------------------------------------------------------
setGeneric("has_cov", function(object) standardGeneric("has_cov"))

#' @describeIn scMethrix-class Checks assay named "counts" (representing coverage) exists
#' @import SummarizedExperiment
setMethod(f = "has_cov", signature = "scMethrix", definition = function(object)   {
  return("counts" %in% SummarizedExperiment::assayNames(object))
})

#---- sampleNames ------------------------------------------------------------------------------------------------------
setGeneric("sampleNames", function(object) standardGeneric("sampleNames"))

#' @describeIn scMethrix-class Gets the sample names stored as rows in [colData()]
setMethod(f = "sampleNames", signature = "scMethrix", definition = function(object)   {
  row.names(colData(object))})

#---- featureNames -----------------------------------------------------------------------------------------------------
setGeneric("featureNames", function(object) standardGeneric("featureNames"))

#' @describeIn scMethrix-class Gets the sample names stored as rows in [rowData()]
setMethod(f = "featureNames", signature = "scMethrix", definition = function(object)   {
  row.names(rowData(object))}
)

#---- Coercion --------------------------------------------------------------------------------------------------------
#' @importClassesFrom minfi GenomicRatioSet
#' @aliases coerce, GenomicRatioSet,scMethrix
#' @exportMethod coerce
setAs("GenomicRatioSet", "scMethrix", function(from) {
  
  if (!is.null(colData(from))) {
    colData <- colData(from)
  } else {
    colData <- data.frame(row.names = colnames(minfi::getBeta(from)))
  }
  
  beta <- minfi::getBeta(from)
  ord <- match(colnames(beta),row.names(colData)) #Ensure colData order consistency
  colData <- colData[ord,,drop=FALSE]
  
  scMethrix(assays = list(score = beta), 
            colData = colData(from), 
            rowRanges = rowRanges(from), 
            is_h5 = FALSE,
            genome = minfi::annotation(from)[["annotation"]])
})

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @exportMethod coerce
setAs("SingleCellExperiment", "scMethrix", function(from) {
  new("scMethrix", as(from, "SingleCellExperiment"))
})

setGeneric("cd", function(object) standardGeneric("cd"))

#---- Metadata functions --------------------------------------------------------------------------
#---- cd ---------------------------------------------------------------------------------------------------------------
#' @describeIn scMethrix-class Same as colData(), but shorter syntax, and will output row names if there is no columns
setMethod(f = "cd", signature = "scMethrix", definition = function(object)   {
  d <- colData(object)
  
  if (ncol(d) == 0) {
    cat(paste("DataFrame with",nrow(d),"rows and 0 columns\n"))
    if (nrow(object) > 10) {
      invisible(sapply(row.names(d)[1:5],function(row) cat(row,"\n")))
      cat("...\n")
      invisible(sapply(row.names(d)[(nrow(d)-5):nrow(d)],function(row) cat(row,"\n")))
    } else {
      invisible(sapply(row.names(d),function(row) cat(row,"\n")))
    }
  } else {
    d
  }
})

#---- rd ---------------------------------------------------------------------------------------------------------------
setGeneric("rd", function(object) standardGeneric("rd"))

#' @describeIn scMethrix-class Same as rowData(), but shorter syntax, and will output row names if there is no columns
setMethod(f = "rd", signature = "scMethrix", definition = function(object)   {
  
  d <- rowData(object)
  
  if (ncol(d) == 0) {
    cat(paste("DataFrame with",nrow(d),"rows and 0 columns\n"))
    if (nrow(object) > 10) {
      invisible(sapply(row.names(d)[1:5],function(row) cat(row,"\n")))
      cat("...\n")
      invisible(sapply(row.names(d)[(nrow(d)-5):nrow(d)],function(row) cat(row,"\n")))
    } else {
      invisible(sapply(row.names(d),function(row) cat(row,"\n")))
    }
  } else {
    d
  }
})

#---- md ---------------------------------------------------------------------------------------------------------------
setGeneric("md", function(object) standardGeneric("md"))

#' @describeIn scMethrix-class Same as metadata(), but shorter syntax
setMethod(f = "md", signature = "scMethrix", definition = function(object)   {
  S4Vectors::metadata(object)
})

#---- is_h5 ------------------------------------------------------------------------------------------------------------
setGeneric("is_h5", function(object) standardGeneric("is_h5"))

#' @describeIn scMethrix-class Checks if \code{\link{scMethrix}} object is an HDF5 object
#' @section is_h5():
#'  This checks the metadata whether the experiment is in HDF5 format. This is found from the metadata attribute \code{is_h5}. A secondary check of all assays also occurs to ensure they are all of appropriate type. An error will be thrown if any assays are the wrong type. This should not occur during normal usage of the package, but may be caused by manual manipulation of assays. If this does occur, [convert_scMethrix()] should be used to restore consistency.
setMethod(f = "is_h5", signature = "scMethrix", definition = function(object)   {
  return (object@metadata$is_h5)
})

#---- .validDims ---------------------------------------------------------------------------
#' Determines if a [scMethrix()] object has valid assay dimensions
#' @param object An [scMethrix()] object
#' @noRd
.validDims <- function(object) {
  # Check the dimensions of assays are correct
  assay_dim <- lapply(assays(object),dim)
  exp_dim <- c(nrow(object),ncol(object))
  match_dim <- sapply(assay_dim, identical, y=exp_dim)
  if (!all(match_dim)) {
    assay_names <- names(which(!match_dim))
    return (paste0("   Wrong dim: Assay",  ifelse(length(assay_names > 1),"s","")," '", paste0(assay_names,collapse="', '"),
                               "' should have dimensions of (",paste0(exp_dim,collapse=", "),")."))
  }
  return (NULL)
}

#---- .validH5 ---------------------------------------------------------------------------
#' Make sure \code{is_h5} variable is present
#' @param object An [scMethrix()] object
#' @noRd
.validH5 <- function(object) {
  if (!"is_h5" %in% names(S4Vectors::metadata(object))) {
    return (paste0("   Missing variable: 'is_h5' is not present in metadata(). Add manually with:\n",
                               "        metadata(object) <- append(metadata(object),list(is_h5 = [TRUE/FALSE]))"))
  }
  return (NULL)
}

#---- .validSamples ---------------------------------------------------------------------------
#' Check if all the sample names are consistent
#' @param object An [scMethrix()] object
#' @noRd
.validSamples <- function(object) {
  assay_cols <- lapply(assays(object),colnames)
  match_cols <- sapply(assay_cols, identical, y=sampleNames(object))
  if (!all(match_cols)) {
    col_names <- names(which(!match_cols))
    return (paste0("   Wrong colNames: Assay",  ifelse(length(col_names > 1),"s","")," '", paste0(col_names,collapse="', '"),
                               "' contain column names that do not match colData()."))
  }
}

#---- .validAssays ---------------------------------------------------------------------------
# Check if all assays are either matrix or HDF5matrix-related types
#' @param object An [scMethrix()] object
#' @noRd
.validAssays <- function(object) {
  exp_type = if (is_h5(object)) c("HDF5Matrix","DelayedMatrix") else "matrix"
  assay_class <- lapply(assays(object),class)
  match_class <- lapply(assay_class, function (a) exp_type %in% a)
  match_class <- sapply(match_class,any)
  if (!all(match_class)) {
    assay_names <- names(which(!match_class))
    return (paste0("   Wrong type: Assay",  ifelse(length(assay_names > 1),"s","")," '", paste0(assay_names,collapse="', '"),
                               "' are of wrong type. Should be ",paste0(exp_type,collapse=" or "),". Recommend using convert_scMethrix() to force type of all assays."))
  }  
  
  return (NULL)
}

#---- .validRedDim ---------------------------------------------------------------------------
# Check if all reduced dim names are consistent
#' @param object An [scMethrix()] object
#' @noRd
.validRedDim <- function(object) {  
  red_rows <- lapply(reducedDims(object),rownames)
  if (!is.empty(red_rows)) {
    match_red <- sapply(red_rows, identical, y=sampleNames(object))
    if (!all(match_red)) {
      red_names <- names(which(!match_red))
      return (paste0("   Wrong rowNames: Reduced dim",  ifelse(length(red_names > 1),"s","")," '", paste0(red_names,collapse="', '"),
                                 "' contain row names that do not match colData()."))
    }
  }
  return(NULL)
}

#---- .validscMethrix ---------------------------------------------------------------------------
#' Determines if a [scMethrix()] object is valid.
#' 
#' Checks for:
#' * 'is_h5' is present in \code{metadata()}
#' * assay dimensions  are consistent with \code{rowData()} and \code{colData()}
#' * sample names in assays (\code{colnames())} are consistent with \code{colData()}
#' * sample names in reduced dims (\code{reducedDims(rownames())}) are consistent with \code{colData()}
#' * assay classes are consistant with \code{is_h5}
#' @param object A [scMethrix()] object
#' @noRd
.validscMethrix <- function(object) {

  errors <- c(
    .validH5(object),
    .validDims(object),
    .validSamples(object),
    .validAssays(object)
  )

  # There is an internal check inside reducedDims in SingleCellExperiment that causes an infinite loop
  #      see: https://github.com/drisso/SingleCellExperiment/blob/master/R/reducedDims.R#L142)
  #      TODO: Capture these errors properly
  # .validRedDim(object)

  if (!is.null(errors)) return (errors)
  
  return (TRUE)
}

S4Vectors::setValidity2("scMethrix", .validscMethrix)

#---- generic_scMethrix_function -----------------------------------------------------------------------------
#' Function used only for inheritance for Roxygen2 documentation. Lists the common function inputs used in the package
#' @param scm [scMethrix()]; a single cell methylation experiment object
#' @param assay string; name of an existing assay. Default = "score"
#' @param new_assay string; name for transformed assay. Default = "new_assay"
#' @param trans closure; The transformation function. Default = mean
#' @param verbose boolean; Flag for outputting function status messages. Default = TRUE 
#' @param n_chunks integer; Number of chunks to split the \code{\link{scMethrix}} object in case it is very large. Default = 1
#' @param n_threads integer; Maximum number of parallel instances. Default = 1
#' @param batch_size integer; The maximum number of elements to process at once.
#' @param h5_dir string; The directory to store HDF5 files. Will be created if it does not exist. Default = NULL
#' @param replace boolean; flag for whether to delete the contents of \code{h5_dir} before saving 
#' @param overlap_type defines the type of the overlap of the CpG sites with the target region. Default value is `within`. For detailed description, see the \code{findOverlaps} function of the \code{\link{IRanges}} package.
#' @param na.rm boolean; flag to remove NA values
generic_scMethrix_function <- function(scm, assay, new_assay, trans, verbose, n_chunks, n_threads, h5_dir, overlap_type, batch_size, replace, na.rm) {}