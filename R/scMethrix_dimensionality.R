#--------------------------------------------------------------------------------------------------------------------------
#' Reduces a matrix
#' @param scm Input \code{\link{scMethrix}} object
#' @param assay The assay to use. Default is 'score'
#' @param top_var Number of variable CpGs to use. Default 1000 Set it to NULL to use all CpGs (which is not recommended due to memory requirements). This option is mutually exclusive with \code{ranges}.
#' @param var Choose between random CpG sites ('rand') or most variable CpGs ('top').
#' @param verbose flag to output messages or not
#' @return Reduced matrix
#' @importFrom stats complete.cases var
#' @examples
#' data('scMethrix_data')
#' reduce_matrix(scMethrix_data)
#' @export
reduce_matrix <- function(scm, assay = "score", var = "top", top_var = 1000, verbose = FALSE) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }
  
  if (!is.numeric(top_var)){
    stop("Parameters top_var must be numeric.")
  }
  
  var_select <- match.arg(var, c("top", "rand"))
  
  if (verbose) message("Generating reduced dataset...")
  
  if (is.null(top_var)) {
    message("All CpGs in the dataset will be used for the reduction")
    meth_sub <- get_matrix(scm = scm, assay = assay, add_loci = FALSE)
  } else {
    
    top_var <- as.integer(as.character(top_var))
    
    if (var_select == "rand") {
      message("Random CpGs within provided GRanges will be used for the reduction")
      ids <- sample(x = seq_along(meth_sub), replace = FALSE, size = min(top_var, nrow(meth_sub)))
      meth_sub <- get_matrix(scm = scm[ids, ], assay = assay, add_loci = FALSE)
    } else {
      message("Taking top ",top_var," most variable CpGs for the reduction")
      meth_sub <- get_matrix(scm = scm, assay = assay, add_loci = FALSE)
      
      if (is_h5(scm)) {
        sds <- DelayedMatrixStats::rowSds(meth_sub, na.rm = TRUE)
      } else {
        sds <- matrixStats::rowSds(meth_sub, na.rm = TRUE)
      }
      
      meth_sub <- meth_sub[order(sds, decreasing = TRUE)[seq_len(min(top_var, nrow(meth_sub)))], ]
    }
  }
  
  # Remove NA
  if (is_h5(scm)) {
    meth_sub <- meth_sub[!DelayedMatrixStats::rowAnyMissings(meth_sub), , drop = FALSE]
  } else {
    meth_sub <- meth_sub[stats::complete.cases(meth_sub), , drop = FALSE]
  }
  
  if (nrow(meth_sub) == 0) {
    stop("Zero loci available post NA removal :(")
  }
  
  return (meth_sub)
}

#--------------------------------------------------------------------------------------------------------------
#' Principal Component Analysis
#' @details Do PCA stuff
#' @inheritParams reduce_matrix 
#' @inheritParams stats::prcomp
#' @return \code{\link{scMethrix}} object with reducedDim 'PCA'
#' @importFrom stats prcomp
#' @seealso [plot_pca()] for plotting
#' @examples
#' data('scMethrix_data')
#' pca_scMethrix(scMethrix_data)
#' @export
#'
pca_scMethrix <- function(scm, assay="score", var = "top", top_var = 1000, verbose = FALSE, ...) {
  
  if (verbose) message("Starting PCA ",start_time())
  
  meth <- reduce_matrix(scm,assay = assay, ...)
  meth <- prcomp(x = t(meth), retx = TRUE, ...)
  
  # Variance explained by PC's
  pc_vars <- meth$sdev^2/sum(meth$sdev^2)
  names(pc_vars) <- colnames(meth$x)
  pc_vars <- round(pc_vars, digits = 2)
  
  reducedDim(scm, "PCA") <- meth$x
  scm@metadata$PCA_vars <- pc_vars

  message("PCA vars (saved in metadata$PCA_vars):")
  print(pc_vars, quote=FALSE, print.gap=2)
  
  if (verbose) message("PCA finished in ",stop_time())

  gc(verbose = FALSE)
    
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Generates UMAP for scMethrix
#' @details Does UMAP stuff
#' @param n_neighbors integer; number of nearest neighbors
#' @inheritParams reduce_matrix 
#' @inheritParams umap::umap
#' @return \code{\link{scMethrix}} object with reducedDim 'UMAP'
#' @import umap
#' @seealso [plot_umap()] for plotting
#' @examples
#' data('scMethrix_data')
#' umap_scMethrix(scMethrix_data)
#' @export
umap_scMethrix <- function(scm, assay="score", var = "top", top_var = 1000, verbose = FALSE, n_neighbors = 15, ...) {
  
  if (verbose) message("Starting UMAP",start_time())
  
  meth <- reduce_matrix(scm,assay = assay, ...)
  umap <- umap(as.matrix(t(meth)),n_neighbors=min(n_neighbors,ncol(scm)), ...)
  
  reducedDim(scm, "UMAP") <- umap$layout
  
  if (verbose) message("UMAP generated in ",stop_time())
  
  gc(verbose = FALSE)
  
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Generates tSNE for scMethrix
#' @details Does tSNE stuff
#' @inheritParams reduce_matrix 
#' @inheritParams Rtsne::Rtsne
#' @return \code{\link{scMethrix}} object with reducedDim 'tSNE'
#' @import Rtsne
#' @seealso [plot_tsne()] for plotting
#' @examples
#' data('scMethrix_data')
#' umap_scMethrix(scMethrix_data)
#' @export
tsne_scMethrix <- function(scm, assay="score", var = "top", top_var = 1000, perplexity = 30, verbose = FALSE, ...) {
  
  if (verbose) message("Starting tSNE",start_time())
  
  meth <- reduce_matrix(scm,assay = assay, ...)
  meth_sub <- Rtsne(as.matrix(t(meth)), perplexity = min(perplexity,floor(ncol(meth)/3)), ...)
  
  reducedDim(scm, "tSNE") <- meth_sub$Y
  
  if (verbose) message("tSNE generated in ",stop_time())
  
  gc(verbose = FALSE)
  
  return(scm)
}