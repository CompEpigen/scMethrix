#--- reduce_scMethrix ----------------------------------------------------------------------------------------
#' Reduces a assay to a representative matrix
#' @details For the purposes of dimensionality reduction, this function selects either random CpGs or those with the highest variability. 
#' @inheritParams generic_scMethrix_function
#' @param n_cpg integer; Number of variable CpGs to use. Default 1000.
#' @param var strning; Choose between random CpG sites ('rand') or most variable CpGs ('top'). Default 'top'
#' @param na.rm boolean; flag to remove NA values
#' @return matrix; the reduced form of the input assay
#' @importFrom stats complete.cases var
#' @examples
#' data('scMethrix_data')
#' nrow(scMethrix_data)
#' scMethrix_data <- reduce_scMethrix(scMethrix_data)
#' nrow(scMethrix_data)
#' @export
reduce_scMethrix <- function(scm, assay = "score", var = c("top", "rand"), n_cpg = 1000, na.rm = FALSE, verbose = FALSE) {

  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  var <- .validateArg(var,reduce_scMethrix)
  .validateType(n_cpg,"integer")
  .validateType(na.rm,"boolean")
  .validateType(verbose,"boolean")
  
  #- Function code -----------------------------------------------------------------------------  
  if (verbose) message("Generating reduced dataset...",start_time())
  
  meth <- get_matrix(scm = scm, assay = assay)
  
  if (var == "rand") {
    if (verbose) message("Selecting ,",n_cpg," random CpGs...")
    ids <- sample(x = 1:nrow(meth), replace = FALSE, size = min(n_cpg, nrow(meth)))
  } else {
    if (verbose) message("Selecting top ",n_cpg," most variable CpGs...")
    sds <- DelayedMatrixStats::rowSds(meth, na.rm = TRUE)
    ids <- head(order(sds,decreasing = TRUE),n_cpg)
  }
  
  ids <- sort(ids)

  if (na.rm) {
    n <- length(ids)
    if (is_h5(scm)) {
      ids <- ids[!DelayedMatrixStats::rowAnyMissings(meth[ids,])]
    } else {
      ids <- ids[stats::complete.cases(meth[ids,])]
    }
  }
  
  if (length(ids) == 0) {
    stop("Zero loci available post NA removal :(")
  } else if (na.rm) {
    message("Removed ",n - length(ids)," CpGs due to NA")
  }

  return (scm[ids,])
}

#--- dim_red_scMethrix --------------------------------------------------------------------------------------
#' Reduces dimensionality (tSNE, UMAP, PCA, or custom)
#' @details Does reduction stuff
#' @param n_components integer; Number of components to use
#' @param n_neighbors integer; number of nearest neighbors for UMAP
#' @param type string; the type of imputation "tSNE","UMAP", or "PCA"
#' @inheritParams generic_scMethrix_function
#' @inheritParams Rtsne::Rtsne
#' @inheritParams umap::umap
#' @inheritParams stats::prcomp
#' @return \code{\link{scMethrix}} object with reducedDim assay
#' @import Rtsne
#' @import umap
#' @importFrom stats prcomp
#' @seealso [plot_dim_red()] for plotting
#' @examples
#' data('scMethrix_data')
#' 
#' @export
dim_red_scMethrix <- function(scm, assay="score", type=c("tSNE","UMAP","PCA"), perplexity = 30, verbose = FALSE, n_components = 2, n_neighbors = 15, ...) {

  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  type <- .validateArg(type,dim_red_scMethrix)
  .validateType(perplexity,"integer")
  .validateType(n_components,"integer")
  .validateType(n_neighbors,"integer")
  .validateType(verbose,"boolean")
  
  #- Function code -----------------------------------------------------------------------------
  if (verbose) message("Starting dimensionality reduction...",start_time())

  meth <- get_matrix(scm,assay = assay)
  
  if (type == "tSNE") {
    
    meth <- Rtsne::Rtsne(as.matrix(t(meth)), perplexity = min(perplexity,floor(ncol(meth)/3)), k = n_components, check_duplicates=F)#, ...)
    
    SingleCellExperiment::reducedDim(scm, "tSNE") <- meth$Y
    
    if (verbose) message("tSNE generated in ",stop_time())
    
  } else if (type == "UMAP") {
    
    umap <- umap::umap(as.matrix(t(meth)),n_neighbors=min(n_neighbors,ncol(scm)),n_components=n_components)#, ...)
    
    SingleCellExperiment::reducedDim(scm, "UMAP") <- umap$layout
    
  } else if (type == "PCA") {
    
    meth <- stats::prcomp(x = as.matrix(t(meth)), retx = TRUE)#, ...)
    
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
  
  gc(verbose = FALSE)
  
  return(scm)
}

