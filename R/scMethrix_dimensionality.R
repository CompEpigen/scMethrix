#-------------------------------------------------------------------------------------------------------------
#' Reduces a assay to a representative matrix
#' @details For the purposes of dimensionality reduction, this function selects either random CpGs or those with the highest variability. 
#' @param scm scMethrix; Input \code{\link{scMethrix}} object
#' @param assay string; The assay to use. Default is 'score'
#' @param top_var integer; Number of variable CpGs to use. Default 1000 Set it to NULL to use all CpGs (which is not recommended due to memory requirements).
#' @param var strning; Choose between random CpG sites ('rand') or most variable CpGs ('top'). Default 'top'
#' @param verbose boolean; flag to output messages or not
#' @param na.rm boolean; flag to remove NA values
#' @return matrix; the reduced form of the input assay
#' @importFrom stats complete.cases var
#' @examples
#' data('scMethrix_data')
#' reduce_cpgs(scMethrix_data)
#' @export
reduce_cpgs <- function(scm, assay = "score", var = "top", top_var = 1000, na.rm = FALSE, verbose = FALSE) {

  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
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
      meth_sub <- get_matrix(scm = scm, assay = assay, add_loci = FALSE)
      ids <- sample(x = 1:nrow(meth_sub), replace = FALSE, size = min(top_var, nrow(meth_sub)))
      meth_sub <- meth_sub[ids, ]
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
  
  if (na.rm) {
      count <- nrow(meth_sub)
    if (is_h5(scm)) {
      meth_sub <- meth_sub[!DelayedMatrixStats::rowAnyMissings(meth_sub), , drop = FALSE]
    } else {
      meth_sub <- meth_sub[stats::complete.cases(meth_sub), , drop = FALSE]
    }
    
    message("Removed ",count - nrow(meth_sub)," CpGs due to NA")
  }
  
  if (nrow(meth_sub) == 0) {
    stop("Zero loci available post NA removal :(")
  }
  
  return (meth_sub)
}

#------------------------------------------------------------------------------------------------------------
#' Reduces dimensionality (tSNE, UMAP, PCA, or custom)
#' @details Does reduction stuff
#' @param n_components integer; Number of components to use
#' @param n_neighbors integer; number of nearest neighbors for UMAP
#' @param type string; the type of imputation "tSNE","UMAP", or "PCA"
#' @inheritParams reduce_cpgs 
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
dim_red_scMethrix <- function(scm, assay="score", type=c("tSNE","UMAP","PCA"), var = "top", top_var = 1000, perplexity = 30, verbose = FALSE, n_components = 2, n_neighbors = 15, ...) {
  
  if (verbose) message("Starting imputation...",start_time())
  
  meth <- reduce_cpgs(scm,assay = assay, var = var, top_var = top_var, verbose = verbose, na.rm = TRUE)
  
  if (type == "tSNE") {
    
    meth_sub <- Rtsne(as.matrix(t(meth)), perplexity = min(perplexity,floor(ncol(meth)/3)), k = n_components)#, ...)
    
    reducedDim(scm, "tSNE") <- meth_sub$Y
    
    if (verbose) message("tSNE generated in ",stop_time())
    
  } else if (type == "UMAP") {
    
    umap <- umap(as.matrix(t(meth)),n_neighbors=min(n_neighbors,ncol(scm)),n_components=n_components)#, ...)
    
    reducedDim(scm, "UMAP") <- umap$layout
    
  } else if (type == "PCA") {
    
    meth <- prcomp(x = as.matrix(t(meth)), retx = TRUE)#, ...)
    
    # Variance explained by PC's
    pc_vars <- meth$sdev^2/sum(meth$sdev^2)
    names(pc_vars) <- colnames(meth$x)
    pc_vars <- round(pc_vars, digits = 2)
    
    reducedDim(scm, "PCA") <- meth$x[,1:n_components]
    scm@metadata$PCA_vars <- pc_vars[1:n_components]
    
    if (verbose) {
      message("PCA vars (saved in metadata$PCA_vars):")
      message(paste(names(pc_vars), pc_vars, sep = ":", collapse = ", "))
      message("PCA finished in ",stop_time())
    }
    
  } else {
    stop("Invalid imputation type specified")
  }
  
  gc(verbose = FALSE)
  
  return(scm)
}

