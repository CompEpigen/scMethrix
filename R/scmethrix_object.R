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

# Create scMethrix obj
#' @param assays 
#'
#' @param colData 
#' @param rowRanges 
#' @param is_hdf5 
#' @param genome_name 
#' @param chrom_size 
#' @param desc 
#' @param h5_dir 
#' @param replace 
#' @param verbose 
#'
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
#' @details 
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
  (x); SummarizedExperiment::assay(x, i="score")}
)

setMethod(f = "sampleNames", signature = "scMethrix", definition = function(object)   {
  row.names(colData(object))}
)

utils::globalVariables(c("sampleNames")) #TODO: find out why this is necessary

# lociNames <- function(x) {
#           row.names(rowData(x))
# }

# setMethod(f = "featureNames", signature = "scMethrix", definition = function(object)   {
#   (object); row.names(colData(object))}
# )

#--- generic_scMethrix_function -----------------------------------------------------------------------------
#' Function used only for inheritence for Roxygen2 documentation. Lists the common function inputs used in the package
#' @param scm \code{\link{scMethrix}}; the single cell methylation experiment
#' @param assay string; name of an existing assay. Default = "score"
#' @param new_assay string; name for transformed assay. Default = "new_assay"
#' @param trans closure; The transformation function. Default = mean
#' @param verbose boolean; Flag for outputting function status messages. Default = TRUE 
#' @param n_chunks integer; Number of chunks to split the \code{\link{scMethrix}} object in case it is very large. Default = 1
#' @param n_threads integer; Maximum number of parallel instances. Default = 1
#' @param batch_size integer; The maximum number of elements to process at once.
#' @param h5_dir string; The directory to use. Will be created if it does not exist. Default = NULL
#' @param replace boolean; flag for whether to delete the contents of h5_dir before saving 
#' @param overlap_type defines the type of the overlap of the CpG sites with the target region. Default value is `within`. For detailed description, see the \code{findOverlaps} function of the \code{\link{IRanges}} package.
#' @param na.rm boolean; flag to remove NA values
generic_scMethrix_function <- function(scm, assay, new_assay, trans, verbose, n_chunks, n_threads, h5_dir, overlap_type, batch_size, replace, na.rm) {}