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
  cat(paste0("An object of class ", class(object), "\n"))
  cat(paste0("   n_CpGs: ", format(nrow(object), big.mark = ","), "\n"))
  cat(paste0("   n_samples: ", ncol(object), "\n"))
  cat(paste0("   assays: ", (paste(SummarizedExperiment::assayNames(object),collapse=", ")),"\n"))
  cat(paste0("   reduced dims: ", (paste(SingleCellExperiment::reducedDimNames(object),collapse=", ")),"\n"))
  cat(paste0("   is_h5: ", is_h5(object), "\n"))
  cat(paste0("   Reference: ", object@metadata$genome, "\n"))
  cat(paste0("   Physical size: ", format(utils::object.size(object), units = "auto"), "\n"))
})

# Create scMethrix obj
create_scMethrix <- function(assays = NULL, colData = NULL, rowRanges = NULL, is_hdf5 = FALSE, 
                             genome_name = "hg19", chrom_size = NULL, desc = NULL, h5_dir = NULL, 
                             replace = FALSE, verbose=TRUE) {

  if (is_hdf5) {

    sse <- SingleCellExperiment::SingleCellExperiment(assays = lapply(assays,function(x) as(x, "HDF5Array")), 
                                                      colData = colData, 
                                                        rowRanges = rowRanges,
                                                        metadata = list(genome = genome_name,
                                                                        chrom_size = chrom_size,
                                                                        descriptive_stats = desc,
                                                                        is_h5 = TRUE))

      if (!is.null(h5_dir)) {
        tryCatch(save_HDF5_scMethrix(scm = sse, h5_dir = h5_dir, replace = replace, verbose = verbose),
                 error = function(e) message(e,"\nThe dataset is not 
                                                         saved. Please save manually using the 
                                                         HDF5Array::saveSummarizedExperiment command."))
      }
      
    } else {
      
      sse <- SingleCellExperiment::SingleCellExperiment(assays = lapply(assays,function(x) as(x, "matrix")),
                                                      colData = colData,
                                                      rowRanges = rowRanges,
                                                      metadata = list(genome = genome_name,
                                                                      chrom_size = chrom_size,
                                                                      descriptive_stats = desc,
                                                                      is_h5 = FALSE))
    }

    return(scMethrix(sse))
}

setMethod(f = "score", signature = "scMethrix", definition = function(x)   {
          (x); SummarizedExperiment::assay(x, i="score")}
)

# setMethod(f = "featureNames", signature = "scMethrix", definition = function(object)   {
#   (object); row.names(colData(object))}
# )

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
generic_scMethrix_function <- function(scm, assay, new_assay, trans, verbose, n_chunks, n_threads, h5_dir, overlap_type) {}