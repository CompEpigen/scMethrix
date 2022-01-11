#' Class scMethrix
#' @description S4 class scMethrix
#' @slot assays A list of two matrices containing 'Methylation' and 'Coverage' information
#' #slot bins A list of matricies for different size of binning for methylation data
#' @slot elementMetadata A DataFrame describing rows in correspoding assay matrices.
#' @slot colData genome: the name of the BSgenome that was used to extract CpGs, isHDF5: is it stored in HDF5 Array format
#' @slot metadata a list of meta data associated with the assays
#' @slot rowRanges A \code{\link{GRanges}} object of the genomic coordinates of CpG sites
#' @slot NAMES NULL
#' @slot int_colData NULL 
#' @slot int_metadata NULL
#' @slot int_elementMetadata NULL 
#' @exportClass scMethrix
#' @importFrom stats median quantile sd
#' @importFrom utils data head write.table menu
#' @importFrom methods is as new
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
scMethrix <- setClass(Class = "scMethrix", contains = "SingleCellExperiment")

setMethod(f = "show", signature = "scMethrix", definition = function(object) {
  
  feature.names <- NULL#row.names(rowData(object))
  # if (!is.null(feature.names)) feature.names <- paste0(substr(feature.names[1:5],1,8),"...",collapse=" l ")
  sample.names <- NULL#row.names(colData(object))
  # if (!is.null(sample.names)) {
  #   sample.names <- substr(row.names(colData(object))[1:5],1,8)
  #   sample.names <- gsub('[^a-zA-Z]*$','',sample.names)
  #   sample.names <- paste0(sample.names,"...",collapse=" l ")
  #   sample.names <- paste0(" (",sample.names,")")
  # }
  
  h5 <- is_h5(object)
  if (h5) h5 = path(score(object))
  
  cat(paste0("An object of class ", class(object), "\n"))
  cat(paste0("   CpGs: ", format(nrow(object), big.mark = ","),feature.names,"\n"))
  cat(paste0("   Samples: ", ncol(object),sample.names,"\n"))
  cat(paste0("   Assays: ", (paste(SummarizedExperiment::assayNames(object),collapse=", ")),"\n"))
  cat(paste0("   Reduced dims: ", (paste(SingleCellExperiment::reducedDimNames(object),collapse=", ")),"\n"))
  cat(paste0("   HDF5: ", h5, "\n"))
  cat(paste0("   Ref.Genome: ", object@metadata$genome, "\n"))
  cat(paste0("   Physical size: ", format(utils::object.size(object), units = "auto"), "\n"))
})

#' Create an scMethrix object
#' @inheritParams generic_scMethrix_function
#' @param assays list of matrices; The assays to include in the experiment
#' @param colData data.frame; the data corresponding to each sample in the assays
#' @param rowRanges GenomicRanges; the genomic loci corresponding to the assays
#' @param is_hdf5 boolean; flag for HDF5 format
#' @param genome_name string; the name of the genome build used
#' @param chrom_size named list of integers; the size of each chromosome in the genomic ranges
#' @param desc named list of strings; list of relevent experiment data
#' @export
create_scMethrix <- function(assays = NULL, colData = NULL, rowRanges = NULL, is_hdf5 = FALSE, 
                             genome_name = "hg19", chrom_size = NULL, desc = NULL, h5_dir = NULL, 
                             replace = FALSE, verbose=TRUE) {
  
  if (is_hdf5) {
    scm <- scMethrix(SingleCellExperiment::SingleCellExperiment(assays = lapply(assays,function(x) as(x, "HDF5Array")), 
                                                      colData = colData, 
                                                      rowRanges = rowRanges,
                                                      metadata = list(genome = genome_name,
                                                                      chrom_size = chrom_size,
                                                                      descriptive_stats = desc,
                                                                      is_h5 = TRUE)))
    if (!is.null(h5_dir)) {
      tryCatch(save_scMethrix(scm = scm, dest = h5_dir, replace = replace, verbose = verbose),
               error = function(e) message(e,"\nThe dataset is not saved. Please save manually using the HDF5Array::saveSummarizedExperiment command."))
    }
    
  } else {
    scm <- scMethrix(SingleCellExperiment::SingleCellExperiment(assays = lapply(assays,function(x) as(x, "matrix")),
                                                      colData = colData,
                                                      rowRanges = rowRanges,
                                                      metadata = list(genome = genome_name,
                                                                      chrom_size = chrom_size,
                                                                      descriptive_stats = desc,
                                                                      is_h5 = FALSE)))
  }
  
  return(scm)
}

#--- as.scMethrix.GRset ----------------------------------------------------------------------------------
#' Converts from a minfi::GRset to an scMethrix object
#' @details Used in converting .IDAT files into scMethrix objects
#' @param GRset GRset; the GRset to convert to scMethrix
#' @param colData data.table; information about the samples
#' @param verbose boolean; be chatty
#' @return An \code{\link{scMethrix}} object
#' @import minfi
#' @examples
#' data('scMethrix_data')
#' @export
as.scMethrix.GRset <- function (GRset, colData = NULL, verbose = verbose) {
  
  if (is.null(colData)) {
    colData <- data.frame(row.names = colnames(getBeta(GRset)))
  } else { # Ensure that colData is in same order as assays
    ord <- match(colnames(getBeta(GRset)),row.names(colData)) 
    colData <- colData[ord,,drop=FALSE]
  }
  
  assays = list(score = minfi::getBeta(GRset))
  rowRanges = rowRanges(GRset)
  genome_name = minfi::annotation(GRset)[["annotation"]]
  
  create_scMethrix(assays = assays, colData = colData, rowRanges = rowRanges, genome_name = genome_name, verbose = verbose)
}

setMethod(f = "score", signature = "scMethrix", definition = function(x)   {
  (x); SummarizedExperiment::assay(x, i="score")})

setMethod(f = "sampleNames", signature = "scMethrix", definition = function(object)   {
  row.names(colData(object))})

utils::globalVariables(c("sampleNames")) #TODO: find out why this is necessary

# lociNames <- function(x) {
#           row.names(rowData(x))
# }

# setMethod(f = "featureNames", signature = "scMethrix", definition = function(object)   {
#   (object); row.names(colData(object))}
# )

#--- Metadata functions --------------------------------------------------------------------------
#' Same as colData(scm), but shorter syntax, and will output row names if there is no columns
#' @inheritParams generic_scMethrix_function
#' @export
cd <- function(scm) {
  
  d <- colData(scm)
  
  if (ncol(d) == 0) {
    cat(paste("DataFrame with",nrow(d),"rows and 0 columns\n"))
    if (nrow(scm) > 10) {
      invisible(sapply(row.names(d)[1:5],function(row) cat(row,"\n")))
      cat("...\n")
      invisible(sapply(row.names(d)[(nrow(d)-5):nrow(d)],function(row) cat(row,"\n")))
    } else {
      invisible(sapply(row.names(d),function(row) cat(row,"\n")))
    }
  } else {
    d
  }
}

#' Same as rowData(scm), but shorter syntax, and will output row names if there is no columns
#' @inheritParams generic_scMethrix_function
#' @export
rd <- function(scm) {
  
  d <- rowData(scm)
  
  if (ncol(d) == 0) {
    cat(paste("DataFrame with",nrow(d),"rows and 0 columns\n"))
    if (nrow(scm) > 10) {
      invisible(sapply(row.names(d)[1:5],function(row) cat(row,"\n")))
      cat("...\n")
      invisible(sapply(row.names(d)[(nrow(d)-5):nrow(d)],function(row) cat(row,"\n")))
    } else {
      invisible(sapply(row.names(d),function(row) cat(row,"\n")))
    }
  } else {
    d
  }
}

#' Same as metadata(scm), but shorter syntax
#' @inheritParams generic_scMethrix_function
#' @export
md <- function(scm) {
  S4Vectors::metadata(scm)
}

#--- is_h5 --------------------------------------------------------------------------------------------------
#' Checks if \code{\link{scMethrix}} object is an HDF5 object
#' @details This checks the metadata whether the experiment is in HDF5 format. This is found from the metadata attribute \code{is_h5}. 
#' 
#' A secondary check of all assays also occurs to ensure they are all of appropriate type. An error will be thrown if any assays are the wrong type. This should not occur during normal usage of the package, but may be caused by manual manipulation of assays. If this does occur, [convert_scMethrix()] should be used to restore consistency.
#' @inheritParams generic_scMethrix_function
#' @return boolean; Whether the object is HDF5
#' @examples
#' data('scMethrix_data')
#' is_h5(scMethrix_data)
#' @export
is_h5 = function(scm) {
  .validateExp(scm)
  
  exp_type = if (scm@metadata$is_h5) c("HDF5Matrix","DelayedMatrix") else "matrix"
  for (name in SummarizedExperiment::assayNames(scm)) {
    if (!any(exp_type %in% class(assay(scm,name)))) 
      stop("Error in scMethrix object. The '",name,"' assay is of type '",
           paste0(class(assay(scm,name)),collapse="', "), "' instead of HDF5Matrix.\nRecommend using convert_scMethrix() to force type of all assays.", call. = FALSE)
  }
  
  return (scm@metadata$is_h5)
}

#--- has_cov ------------------------------------------------------------------------------------------------
#' Checks if [scMethrix()] object has a coverage matrix.
#' @details This check for the existence of a \code{counts} matrix in the object
#' @inheritParams generic_scMethrix_function
#' @return boolean; Whether the object has a coverage matrix
#' @import SummarizedExperiment
#' @examples
#' data('scMethrix_data')
#' has_cov(scMethrix_data)
#' @export
has_cov = function(scm) {
  .validateExp(scm)
  return("counts" %in% SummarizedExperiment::assayNames(scm))
}

#--- generic_scMethrix_function -----------------------------------------------------------------------------
#' Function used only for inheritence for Roxygen2 documentation. Lists the common function inputs used in the package
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