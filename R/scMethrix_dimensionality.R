#--- reduce_scMethrix ----------------------------------------------------------------------------------------
#' Reduces a assay to a representative matrix
#' @details For the purposes of dimensionality reduction, this function selects either random CpGs or those with the highest variability. 
#' @inheritParams generic_scMethrix_function
#' @param n_cpg integer; Number of variable CpGs to reduce to. Default 1000.
#' @param var string; Choose between random CpG sites ('rand') or most variable CpGs ('top'). Default 'top'. Seed for random sampling = n_cpg.
#' @param na.rm boolean; flag to remove features with any NA values. Final number of CpGs will likely be less than n_cpg. For 'random', Will iterate up to 10 times trying to find random cpgs.
#' @return matrix; the reduced form of the input assay
#' @importFrom stats complete.cases var
#' @examples
#' data('scMethrix_data')
#' nrow(scMethrix_data)
#' scMethrix_data <- reduce_scMethrix(scMethrix_data)
#' nrow(scMethrix_data)
#' @export
reduce_scMethrix <- function(scm, assay = "score", var = c("top", "random"), n_cpg = 1000, na.rm = FALSE, verbose = FALSE) {

  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  var <- .validateArg(var,reduce_scMethrix)
  .validateType(n_cpg,"integer")
  .validateType(na.rm,"boolean")
  .validateType(verbose,"boolean")
  
  set.seed(n_cpg)
  
  if (verbose) message("Selecting ",n_cpg," ",var," CpGs...",start_time())
  
  #- Function code -----------------------------------------------------------------------------  

  n_cpg <- min(n_cpg,nrow(scm))

  meth <- get_matrix(scm = scm, assay = assay)
  
  if (na.rm)
    ids <- which(!DelayedMatrixStats::rowAnyNAs(meth))
  else {
    ids <- 1:nrow(scm)
  }

  if (length(ids) > 0) {
    if (var == "random") {
      ids <- sample(x = ids, replace = FALSE, size = n_cpg)
    } else {
      sds <- DelayedMatrixStats::rowSds(meth[ids,], na.rm = TRUE)
      ids <- head(order(sds,decreasing = TRUE),n_cpg)
    }
    
  }
  
  if (length(ids) == 0) {
    stop("No CpGs left after reduction. No non-zero containing rows were found.")
  } else if (na.rm) {
    message("Removed ",n_cpg - length(ids)," CpGs due to NA")
  }

  return (scm[sort(ids),])
}

#--- dim_red_scMethrix --------------------------------------------------------------------------------------
#' Reduces dimensionality (tSNE, UMAP, PCA, or custom)
#' @details Does reduction stuff
#' 
#' The dimensionality reduction algorithms used do not allow for NA values, so the input matrix will likely need to be imputed or binned, among others. 
#' 
#' Certain components have hard-coded minimums for the dimensionality reduction to occur (only a concern from very small sample sets, about 30 or so). These will be automatically enforced:
#' * tSNE: perplexity >= floor(ncol(scm)/3)
#' * UMAP: n_neighbors >= ncol(scm)
#' @param plot boolean; Plot after calculating
#' @param n_components integer; Number of components to use
#' @param n_neighbors integer; number of nearest neighbors for UMAP
#' @param type string; the type of imputation "tSNE","UMAP", or "PCA"
#' @param ... additional arguements for any of the imputation functions
#' @inheritParams generic_scMethrix_function
#' @inheritParams Rtsne::Rtsne
#' @inheritParams umap::umap
#' @inheritParams stats::prcomp
#' @return \code{\link{scMethrix}} object with reducedDim assay
#' @import Rtsne
#' @import umap
#' @importFrom stats prcomp
#' @seealso [plot_dim_red()] for plotting, [Rtsne::Rtsne()] for Rtsne, [umap::umap()] for UMAP, [stats::prcomp()] for PCA
#' @examples
#' data('scMethrix_data')
#' #scMethrix_data <- transform_assay(scMethrix_data,assay="score",new_assay="fill",trans = fill)
#' #dim_red_scMethrix(scMethrix_data, assay="fill", type="PCA")
#' @export
dim_red_scMethrix <- function(scm, assay="score", type=c("tSNE","UMAP","PCA"), plot = F, perplexity = 30, verbose = FALSE, n_components = 2, n_neighbors = 15, ...) {

  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  type <- .validateArg(type,dim_red_scMethrix)
  .validateType(perplexity,"integer")
  .validateType(n_components,"integer")
  .validateType(n_neighbors,"integer")
  .validateType(verbose,"boolean")
  
  meth <- get_matrix(scm,assay = assay)
  if (anyNA(get_matrix(scm,assay = assay))) stop("Assay matrix cannot contain NAs. You must impute or otherwise fill these values.",call. = FALSE)
  
  #- Function code -----------------------------------------------------------------------------
  if (verbose) message("Starting dimensionality reduction...",start_time())

  if (type == "tSNE") {
    
    meth <- Rtsne::Rtsne(as.matrix(t(meth)), perplexity = min(perplexity,floor(ncol(meth)/3)), dims = n_components, check_duplicates=F, ...)
    SingleCellExperiment::reducedDim(scm, "tSNE") <- meth$Y
  
  } else if (type == "UMAP") {
    
    umap <- umap::umap(as.matrix(t(meth)),n_neighbors=min(n_neighbors,ncol(scm)),n_components=n_components, ...)
    SingleCellExperiment::reducedDim(scm, "UMAP") <- umap$layout
    
  } else if (type == "PCA") {
    
    meth <- stats::prcomp(x = as.matrix(t(meth)), retx = TRUE, ...)
    
    # Variance explained by PC's
    pc_vars <- meth$sdev^2/sum(meth$sdev^2)
    names(pc_vars) <- colnames(meth$x)
    pc_vars <- round(pc_vars, digits = 2)
    
    SingleCellExperiment::reducedDim(scm, "PCA") <- meth$x[,1:n_components]
    scm@metadata$PCA_vars <- pc_vars[1:n_components]
    
    if (verbose) {
      message("PCA vars (saved in metadata$PCA_vars):")
      message(paste(names(pc_vars), pc_vars, sep = ":", collapse = ", "))
      message("PCA finished in ",stop_time())
    }
    
  }
  
  if (plot) plot_dim_red(scm,type)
  if (verbose) message("Dimensionality reduction finished in ",stop_time())
  gc(verbose = FALSE)
  
  return(scm)
}

