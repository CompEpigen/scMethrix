#--------------------------------------------------------------------------------------------------------------------------
#' Principal Component Analysis
#'
#' @param scm Input \code{\link{scMethrix}} object
#' @param top_var Number of variable CpGs to use. Default 1000 Set it to NULL to use all CpGs (which is not recommended due to memory requirements). This option is mutually exclusive with \code{ranges}.
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param pheno Column name of colData(m). Default NULL. Will be used as a factor to color different groups
#' @param var Choose between random CpG sites ('rand') or most variable CpGs ('top').
#' @param do_plot Should a plot be generated?
#' @param n_pc Default 2.
#' @importFrom stats complete.cases prcomp
#' @importFrom graphics axis legend lines mtext par plot title
#' @return PCA results
#' @examples
#' data('scMethrix_data')
#' pca_scMethrix(scMethrix_data, do_plot = FALSE)
#' @export
#'
pca_scMethrix <- function(scm, var = "top", top_var = 1000, ranges = NULL,
                          pheno = NULL, do_plot = TRUE, n_pc = 2) {
  var_select <- match.arg(var, c("top", "rand"))
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!(is.numeric(top_var) & is.numeric(n_pc))){
    stop("Either top_var or n_pc variables are not numeric.")
  }
  
  ## subset based on the input ranges
  if (!is.null(ranges)) {
    message("GenomicRanges will be used for the PCA")
    meth_sub <- subset_scMethrix(scm = scm, regions = ranges)
    meth_sub <- get_matrix(scm = meth_sub, type = "score", add_loci = FALSE)
  }
  
  if (is.null(top_var)) {
    message("All CpGs in the dataset will be used for the PCA")
    if (is.null(ranges)) {
      meth_sub <- get_matrix(scm = scm, type = "score", add_loci = FALSE)
    }
  } else {
    if (!is.numeric(top_var)){
      stop("top_var must be numeric.")
    }
    top_var <- as.integer(as.character(top_var))
    if (var_select == "rand") {
      if (!is.null(ranges)) {
        message("Random CpGs within provided GRanges will be used for the PCA")
        ids <- sample(x = seq_along(meth_sub), replace = FALSE,
                      size = min(top_var, nrow(meth_sub)))
      } else {
        message("Random CpGs will be used for the PCA")
        ids <- sample(x = seq_along(scm), replace = FALSE, size = as.integer(as.character(min(top_var,
                                                                                              nrow(scm)))))
      }
      meth_sub <- get_matrix(scm = scm[ids, ], type = "score", add_loci = FALSE)
    } else {
      if (!is.null(ranges)) {
        if (is_h5(scm)) {
          sds <- DelayedMatrixStats::rowSds(meth_sub, na.rm = TRUE)
        } else {
          sds <- matrixStats::rowSds(meth_sub, na.rm = TRUE)
        }
        meth_sub <- meth_sub[order(sds, decreasing = TRUE)[seq_len(min(top_var,
                                                                       nrow(meth_sub)))], ]
      } else {
        meth_sub <- get_matrix(scm = scm,#order_by_sd(scm)[seq_len(min(top_var, nrow(scm)))],
                               type = "score", add_loci = FALSE)
      }
    }
  }
  
  # Remove NA
  meth_sub <- meth_sub[complete.cases(meth_sub), , drop = FALSE]
  if (nrow(meth_sub) == 0) {
    stop("Zero loci available post NA removal :(")
  }
  
  meth_pca <- prcomp(x = t(meth_sub), retx = TRUE)
  
  n_pc <- ncol(meth_pca$x)
  
  # Variance explained by PC's
  pc_vars <- meth_pca$sdev^2/sum(meth_pca$sdev^2)
  names(pc_vars) <- colnames(meth_pca$x)
  pc_vars <- round(pc_vars, digits = 2)
  
  #-----------------------------------------------------------------------------------------------------------------------
  # Draw cumulative variance explained by PCs
  
  if (do_plot) {
    par(bty = "n", mgp = c(2.5, 0.5, 0), mar = c(3, 4, 2, 2) + 0.1,
        tcl = -0.25, las = 1)
    plot(pc_vars, type = "h", col = "red", xlab = "", ylab = "variance Explained",
         ylim = c(0, 1), yaxs = "i")
    mtext(side = 1, "Principal component", line = 2)
    cum_var <- cumsum(meth_pca$sdev^2)/sum(meth_pca$sdev^2) * meth_pca$sdev[1]^2/sum(meth_pca$sdev^2)
    lines(cumsum(cum_var), type = "s")
    axis(side = 4, at = pretty(c(0, 1)), labels = pretty(c(0, 1)))
    legend("topright", col = c("red", "black"), lty = 1, c("Per PC", "Cumulative"), bty = "n")
    #lines(x = c(length(meth_pca$sdev), n_pc, n_pc), y = c(cum_var[n_pc], cum_var[n_pc], 0), lty = 3)
    title(main = paste0("Variance explained by ", n_pc, " PC: ", round(sum(c(meth_pca$sdev^2/sum(meth_pca$sdev^2))[seq_len(n_pc)]),
                                                                       digits = 2)), adj = 0)
  }
  #-----------------------------------------------------------------------------------------------------------------------
  
  results <- list(PC_matrix = meth_pca$x, var_explained = pc_vars)
  
  if (do_plot) {
    plot_pca(pca_res = results, scm = scm, col_anno = pheno)
  }
  
  gc(verbose = FALSE)
  
  return(results)
}


#------------------------------------------------------------------------------------------------------------
#' Generates UMAP for scMethrix
#' @details Does UMAP stuff
#' @param scm A \code{\link{scMethrix}} object
#' @return An \code{\link{scMethrix}} object
#' @import umap
#' @import ggplot2
#' @examples
#' data('scMethrix_data')
#' @export
umap_scMethrix <- function(scm) {
  
  # x <- y <- NULL
  # 
  # start_time()
  # 
  # scm.bin <- transform_assay(scm.small,assay="score",name="binary",trans=binarize)
  # 
  # start_time()
  # scm.umap <- umap(as.matrix(get_matrix(scm.bin,type="bin")),n_neighbors=min(100,ncol(scm.bin)))
  # stop_time()
  # beep()
  # 
  # stop_time()
  # 
  # df <- data.frame(x = scm.umap$layout[,1],
  #                  y = scm.umap$layout[,2])
  # 
  # ggplot(df, aes(x, y)) +
  #   geom_point()
}