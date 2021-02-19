#' Class scMethrix
#' @description S4 class scMethrix
#' @slot assays A list of two matrices containing 'Methylation' and 'Coverage' information
#' @slot bins A list of matricies for different size of binning for methylation data
#' @slot elementMetadata A DataFrame describing rows in correspoding assay matrices.
#' @slot colData genome: the name of the BSgenome that was used to extract CpGs, isHDF5: is it stored in HDF5 Array format
#' @slot metadata a list of meta data associated with the assays
#' @slot NAMES NULL
#' @exportClass scMethrix
#' @importFrom graphics axis legend lines mtext par plot title
#' @importFrom stats complete.cases cov density median prcomp quantile sd
#' @importFrom utils browseURL
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
scMethrix <- setClass(Class = "scMethrix", contains = "SingleCellExperiment")

setMethod(f = "show", signature = "scMethrix", definition = function(object) {
  cat(paste0("An object of class ", class(object), "\n"))
  cat(paste0("   n_CpGs: ", format(nrow(object), big.mark = ","), "\n"))
  cat(paste0("   n_samples: ", ncol(object), "\n"))
  cat(paste0("   assays: ", assayNames(object),"\n"))
  cat(paste0("    is_h5: ", is_h5(object), "\n"))
  cat(paste0("   Reference: ", object@metadata$genome, "\n"))
})

# Create scMethrix obj
create_scMethrix <- function(methyl_mat = NULL, colData = NULL, rowRanges = NULL, is_hdf5 = FALSE, genome_name = "hg19",
                           chrom_sizes = NULL, desc = NULL) {

    if (is_hdf5) {

      sse <-






    } else {

      sse <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as(methyl_mat, "sparseMatrix")),
                                                      colData = colData,
                                                      rowRanges = rowRanges,
                                                      metadata = list(genome = genome_name,
                                                                      chrom_sizes = chrom_sizes,
                                                                      descriptive_stats = desc))
    }

    return(scMethrix(rse))
}
