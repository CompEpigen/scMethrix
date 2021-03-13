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
  cat(paste0("   on_disk: ", object@metadata$on_disk, "\n"))
  cat(paste0("   Reference: ", object@metadata$genome, "\n"))
  cat(paste0("   Physical size: ", format(object.size(object), units = "auto"), "\n"))
})

# Create scMethrix obj
create_scMethrix <- function(methyl_mat = NULL, colData = NULL, rowRanges = NULL, is_hdf5 = FALSE, genome_name = "hg19",
                           chrom_sizes = NULL, desc = NULL) {
    if (is_hdf5) {

      sse <- SingleCellExperiment::SingleCellExperiment(assays = list(score = as(methyl_mat, "HDF5Array")),
                                                        colData = colData,
                                                        rowRanges = rowRanges,
                                                        metadata = list(genome = genome_name,
                                                                        chrom_sizes = chrom_sizes,
                                                                        descriptive_stats = desc,
                                                                        is_hdf5 = TRUE))
      
      if (!is.null(h5_dir)) {
        tryCatch(HDF5Array::saveHDF5SummarizedExperiment(x = sse, dir = h5_dir,
                                                         replace = TRUE), error = function(e)
                                                           message("The dataset is not saved. Please save manually, using the HDF5Array::saveSummarizedExperiment command. "))
      }
      
    } else {

      sse <- SingleCellExperiment::SingleCellExperiment(assays = list(score = methyl_mat),
                                                      #colData = colData,
                                                      rowRanges = rowRanges,
                                                      metadata = list(genome = genome_name,
                                                                      chrom_sizes = chrom_sizes,
                                                                      descriptive_stats = desc,
                                                                      is_hdf5 = FALSE))
    }

    return(scMethrix(sse))
}

