#--------------------------------------------------------------------------------------------------------------------------
#' Reduces a matrix
#' @param scm Input \code{\link{scMethrix}} object
#' @param assay The assay to use. Default is 'score'
#' @param top_var Number of variable CpGs to use. Default 1000 Set it to NULL to use all CpGs (which is not recommended due to memory requirements). This option is mutually exclusive with \code{ranges}.
#' @param var Choose between random CpG sites ('rand') or most variable CpGs ('top').
#' @param verbose flag to output messages or not
#' @return Reduced matrix
#' @examples
#' data('scMethrix_data')
#' reduce_matrix(scMethrix_data)
#' @export
reduce_matrix <- function(scm, assay = "score", var = "top", top_var = 1000,
                            pheno = NULL, verbose = FALSE) {

  if (!is.numeric(top_var)){
    stop("Parameters top_var must be numeric.")
  }
    
  var_select <- match.arg(var, c("top", "rand"))
  
  if (verbose) message("Generating reduced dataset...",start_time())
  
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
    meth_sub <- meth_sub[!DelayedMatrixStats ::rowAnyMissings(meth_sub), , drop = FALSE]
  } else {
    meth_sub <- meth_sub[complete.cases(meth_sub), , drop = FALSE]
  }
  
  if (nrow(meth_sub) == 0) {
    stop("Zero loci available post NA removal :(")
  }
  
  return (meth_sub)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Principal Component Analysis
#'
#' @param scm Input \code{\link{scMethrix}} object
#' @param assay The assay to use. Default is 'score'
#' @param top_var Number of variable CpGs to use. Default 1000. Set it to NULL to use all CpGs (which is not recommended due to memory requirements). This option is mutually exclusive with \code{ranges}.
#' @param var Choose between random CpG sites ('rand') or most variable CpGs ('top').
#' @param verbose flag to output messages or not
#' @return \code{\link{scMethrix}} object with reducedDim 'PCA'
#' @examples
#' data('scMethrix_data')
#' pca_scMethrix(scMethrix_data)
#' @export
#'
pca_scMethrix <- function(scm, assay="score", var = "top", top_var = 1000,
                          verbose = FALSE) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }
  
  meth_sub <- reduce_matrix(scm,assay = assay, var = var, top_var = top_var, pheno = pheno, verbose = verbose)
  
  meth_pca <- prcomp(x = t(meth_sub), retx = TRUE)
  
    # Variance explained by PC's
  pc_vars <- meth_pca$sdev^2/sum(meth_pca$sdev^2)
  names(pc_vars) <- colnames(meth_pca$x)
  pc_vars <- round(pc_vars, digits = 2)
    
  reducedDim(scm, "PCA") <- meth_pca$x
  scm@metadata$PCA_vars <- pc_vars
  
  if (verbose) message("PCA finished in ",stop_time())
    
  message("PCA vars (saved in metadata$PCA_vars):")
  print(pc_vars, quote=FALSE, print.gap=2)
    
  gc(verbose = FALSE)
  
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Generates UMAP for scMethrix
#' @details Does UMAP stuff
#' @param scm Input \code{\link{scMethrix}} object
#' @param assay The assay to use. Default is 'score'
#' @param top_var Number of variable CpGs to use. Default 1000. Set it to NULL to use all CpGs (which is not recommended due to memory requirements). This option is mutually exclusive with \code{ranges}.
#' @param var Choose between random CpG sites ('rand') or most variable CpGs ('top').
#' @param verbose flag to output messages or not
#' @return \code{\link{scMethrix}} object with reducedDim 'umap'
#' @import umap
#' @examples
#' data('scMethrix_data')
#' umap_scMethrix(scMethrix_data)
#' @export
umap_scMethrix <- function(scm, assay="score", n_neighbors = 20, var = "top", top_var = 1000,
                           verbose = FALSE) {
  x <- y <- NULL
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }

  meth_sub <- reduce_matrix(scm,assay = assay, var = var, top_var = top_var, pheno = pheno, verbose = verbose)
  umap <- umap(as.matrix(t(meth_sub)),n_neighbors=min(n_neighbors,ncol(scm.bin)))

  reducedDim(scm, "umap") <- umap$layout
  
  return(scm)
}