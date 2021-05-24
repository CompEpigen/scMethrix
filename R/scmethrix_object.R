#' Class scMethrix
#' @description S4 class scMethrix
#' @slot assays A list of two matrices containing 'Methylation' and 'Coverage' information
#' #slot bins A list of matricies for different size of binning for methylation data
#' @slot elementMetadata A DataFrame describing rows in correspoding assay matrices.
#' @slot colData genome: the name of the BSgenome that was used to extract CpGs, isHDF5: is it stored in HDF5 Array format
#' @slot metadata a list of meta data associated with the assays
#' @slot rowRanges A \code{\link{GRanges)} object of the genomic coordinates of CpG sites
#' @slot NAMES NULL
#' @slot int_colData NULL 
#' @slot int_metadata NULL
#' @slot int_elementMetadata NULL 
#' @exportClass scMethrix
#' @importFrom graphics axis legend lines mtext par plot title
#' @importFrom stats complete.cases cov density median prcomp quantile sd
#' @import utils methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
scMethrix <- setClass(Class = "scMethrix", contains = "SingleCellExperiment")

setMethod(f = "show", signature = "scMethrix", definition = function(object) {
  cat(paste0("An object of class ", class(object), "\n"))
  cat(paste0("   n_CpGs: ", format(nrow(object), big.mark = ","), "\n"))
  cat(paste0("   n_samples: ", ncol(object), "\n"))
  cat(paste0("   assays: ", (paste(SummarizedExperiment::assayNames(object),collapse=", ")),"\n"))
  cat(paste0("   is_h5: ", is_h5(object), "\n"))
  cat(paste0("   Reference: ", object@metadata$genome, "\n"))
  cat(paste0("   Physical size: ", format(utils::object.size(object), units = "auto"), "\n"))
})

# Create scMethrix obj
create_scMethrix <- function(assays = NULL, colData = NULL, rowRanges = NULL, is_hdf5 = FALSE, 
                             genome_name = "hg19", chrom_sizes = NULL, desc = NULL, h5_dir = NULL, 
                             replace = FALSE, verbose=TRUE) {
    
  if (is_hdf5) {

    sse <- SingleCellExperiment::SingleCellExperiment(assays = lapply(assays,function(x) as(x, "HDF5Array")), 
                                                      colData = colData, 
                                                        rowRanges = rowRanges,
                                                        metadata = list(genome = genome_name,
                                                                        chrom_sizes = chrom_sizes,
                                                                        descriptive_stats = desc,
                                                                        is_h5 = TRUE, 
                                                                        has_cov = ("coverage" %in% names(assays))))
      #TODO: Cannot save to same directory input files exist in
      if (!is.null(h5_dir)) {
        message("Writing to disk...",start_time())
        
        tryCatch(HDF5Array::saveHDF5SummarizedExperiment(x = sse, dir = h5_dir, replace = replace, 
                                                         chunkdim = c(length(rowRanges),1), verbose=verbose), 
                                                         error = function(e) message(e,"\nThe dataset is not 
                                                         saved. Please save manually using the 
                                                         HDF5Array::saveSummarizedExperiment command."))
        message("Written in ",stop_time())
        }
      
    } else {
      
      sse <- SingleCellExperiment::SingleCellExperiment(assays = lapply(assays,function(x) as(x, "matrix")),
                                                      colData = colData,
                                                      rowRanges = rowRanges,
                                                      metadata = list(genome = genome_name,
                                                                      chrom_sizes = chrom_sizes,
                                                                      descriptive_stats = desc,
                                                                      is_h5 = FALSE, 
                                                                      has_cov = ("coverage" %in% names(assays))))
    }

    return(scMethrix(sse))
}

