#' Format \code{\link{scMethrix}} matrix to long form data for plotting
#' @inheritParams generic_scMethrix_function
#' @param n_cpgs Use these many random CpGs for plotting. Default 25000. Set it to \code{NULL} to use all - which can be memory expensive.
#' @param ranges genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param pheno Row name of colData(m). Will be used as a factor to color different groups
#' @return 'Long' matrix for methylation
#' @export
prepare_plot_data <- function(scm = NULL, ranges = NULL, n_cpgs = 25000, pheno = NULL){
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  if (!is.null(n_cpgs)){
    if (!is.numeric(n_cpgs)){
      stop("n_cpgs must be numeric.")
    }
  }
  
  if (!is.null(ranges)) {
    meth_sub <- subset_scMethrix(scm = scm, regions = ranges)
    if (!is.null(n_cpgs)) {
      message("Randomly selecting ", n_cpgs, " sites")
      ids <- sample(x = seq_along(meth_sub), replace = FALSE, size = min(n_cpgs,
                                                                         nrow(meth_sub)))
      meth_sub <- get_matrix(scm = meth_sub[ids, ], add_loci = FALSE)
    } else {
      meth_sub <- get_matrix(scm = meth_sub, add_loci = FALSE)
    }
  } else if (!is.null(n_cpgs)) {
    message("Randomly selecting ", n_cpgs, " sites")
    
    ids <- sample(x = seq_along(scm), replace = FALSE, size = min(n_cpgs,
                                                                nrow(scm)))
    meth_sub <- get_matrix(scm = scm[ids, ], add_loci = FALSE)
  } else {
    meth_sub <- get_matrix(scm = scm, add_loci = FALSE)
  }
  
  
  if (!is.null(pheno)) {
    if (pheno %in% colnames(colData(scm))) {
      colnames(meth_sub) <- as.character(scm@colData[, pheno])
    } else {
      stop("Please provide a valid phenotype annotation column.")
    }
  }
  
  meth_sub <- as.data.frame(meth_sub)
  data.table::setDT(x = meth_sub)
  plot.data <- suppressWarnings(data.table::melt(meth_sub))
  colnames(plot.data) <- c("variable", "Meth")
  
  gc(verbose = FALSE)
  return(plot.data)
  
}

#' Getter for plot palette colors
#' @param n_row Number of colors
#' @param col_palette String for RColorBrewer palette name  
#' @return RColorBrewer palette
get_palette <- function(n_row, col_palette){
  
  if (n_row == 0) {
    stop("Zero colors present in the palette")
  }
  
  if (!col_palette %in% rownames(RColorBrewer::brewer.pal.info)){
    stop("Please provide a valid RColorBrewer palettte. Possible values are: ", paste0(rownames(RColorBrewer::brewer.pal.info)), sep=", ")
  }
  color_pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[col_palette,
                                                                                                  "maxcolors"], col_palette))(n_row)
  return(color_pal)
}

#' Violin Plot for \eqn{\beta}-Values
#' @inheritParams prepare_plot_data
#' @param col_palette Name of the RColorBrewer palette to use for plotting.
#' @return ggplot2 object
#' @export
#' @import ggplot2
#' @examples
#' data('scMethrix_data')
#' plot_violin(scm = scMethrix_data)
plot_violin <- function(scm = NULL, ranges = NULL, n_cpgs = 25000, pheno = NULL,
                        col_palette = "RdYlGn", show_legend = TRUE) {
  variable <- Meth <- NULL
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (is.null(ranges)) ranges = rowRanges(scm)
  
  plot.data <- prepare_plot_data(scm=scm, ranges = ranges, n_cpgs = n_cpgs, pheno = pheno)
  
  col_palette <- get_palette(ncol(scm), col_palette)
  # generate the violin plot
  p <- ggplot2::ggplot(plot.data, ggplot2::aes(x = variable, y = Meth, fill = variable)) + 
    ggplot2::geom_violin(alpha = 0.8, show.legend = show_legend) + ggplot2::theme_classic(base_size = 14) +
    ggplot2::scale_fill_manual(values = col_palette) +
    ggplot2::xlab(pheno) + ggplot2::ylab(expression(beta * "-Value")) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12,
                                                                     colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),
          axis.title.y = element_blank(), legend.title = element_blank())
  
  return(p + scMethrix_theme())
}

#--------------------------------------------------------------------------------------------------------------------------
#' Density Plot of \eqn{\beta}-Values
#'
#' @inheritParams plot_violin
#' @return ggplot2 object
#' @export
#' @examples
#' data('scMethrix_data')
#' plot_density(scm = scMethrix_data)
plot_density <- function(scm = NULL, ranges = NULL, n_cpgs = 25000, pheno = NULL,
                         col_palette = "RdYlGn", show_legend = TRUE) {
  
  variable <- Meth <- NULL
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (is.null(ranges)) ranges = rowRanges(scm)
  
  plot.data <- prepare_plot_data(scm=scm, ranges = ranges, n_cpgs = n_cpgs, pheno = pheno)
  col_palette <- get_palette(ncol(scm), col_palette)
  
  # generate the density plot
  p <- ggplot2::ggplot(plot.data, ggplot2::aes(Meth, color = variable)) +
    geom_density(lwd = 1, position = "identity", show.legend = show_legend) + ggplot2::theme_classic() +
    ggplot2::xlab("Methylation") + ggplot2::theme_classic(base_size = 14) +
    ggplot2::scale_fill_manual(values = col_palette) +
    ggplot2::xlab(expression(beta * "-Value")) + theme(axis.title.x = element_blank(),
                                                       axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12,
                                                                                                                                           colour = "black"), axis.title.y = element_blank(), legend.title = element_blank())
  
  gc(verbose = FALSE)
  
  return(p + scMethrix_theme())
}

#--------------------------------------------------------------------------------------------------------------------------
#' Coverage QC Plots
#' @inheritParams plot_violin
#' @param perGroup Color the plots in a sample-wise manner?
#' @param lim Maximum coverage value to be plotted.
#' @param type Choose between 'hist' (histogram) or 'dens' (density plot).
#' @param size.lim The maximum number of observarions (sites*samples) to use. If the dataset is larger that this,
#' random sites will be selected from the genome.
#' @return ggplot2 object
#' @examples
#' data('scMethrix_data')
#' plot_coverage(scm = scMethrix_data)
#' @export
plot_coverage <- function(scm = NULL, type = c("hist", "dens"), pheno = NULL, perGroup = FALSE,
                          lim = 100, size.lim = 1e+06, col_palette = "RdYlGn", show_legend = TRUE) {
  
  value <- variable <- NULL
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  colors_palette <- get_palette(ncol(scm), col_palette)
  
  type <- match.arg(arg = type, choices = c("hist", "dens"), several.ok = FALSE)
  
  if (nrow(scm) > size.lim) {
    message("The dataset is bigger than the size limit. A random subset of the object will be used that contains ~",
            size.lim, " observations.")
    n_rows <- trunc(size.lim/nrow(scm@colData))
    sel_rows <- sample(seq_len(nrow(scm@elementMetadata)), size = n_rows,
                       replace = FALSE)
    
    meth_sub <- get_matrix(scm = scm[sel_rows, ], assay = "counts",
                                    add_loci = FALSE)
    
  } else {
    meth_sub <- get_matrix(scm = scm, assay = "counts", add_loci = FALSE)
  }
  
  if (perGroup) {
    if (is.null(pheno)) {
      stop("For group based plotting, provide group information using the pheno argument.")
    }
    if (pheno %in% rownames(colData(scm)) == 0) {
      stop("Phenotype annotation cannot be found in colData(m).")
    }
    colnames(meth_sub) <- scm@colData[, pheno]
  }
  
  meth_sub <- as.data.frame(meth_sub)
  data.table::setDT(x = meth_sub)
  plot.data <- suppressWarnings(expr = data.table::melt(meth_sub))
  
  plot.data <- plot.data[value <= lim, ]
  
  # generate the plots
  if (!perGroup) {
    if (type == "dens") {
      p <- ggplot2::ggplot(plot.data, aes(value, color = variable)) +
        ggplot2::geom_density(alpha = 0.5, adjust = 1.5, lwd = 1, show.legend = show_legend,
                              position = "identity") + ggplot2::theme_classic() + ggplot2::xlab("Coverage") +
        ggplot2::scale_fill_manual(values = colors_palette)
      
    } else if (type == "hist") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(value, fill = variable)) + 
        ggplot2::geom_histogram(alpha = 0.6, binwidth = 1, color = "black", show.legend = show_legend) + 
        ggplot2::theme_classic() +
        ggplot2::xlab("Coverage")+
        ggplot2::scale_fill_manual(values = colors_palette)
      # print(p)
    }
  } else {
    if (type == "dens") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(value, color = variable)) +
        ggplot2::geom_density(alpha = 0.6, adjust = 1.5, lwd = 1, show.legend = show_legend,
                              position = "identity") + ggplot2::theme_classic() + ggplot2::xlab("Coverage") +
        ggplot2::labs(fill = "Groups") +
        ggplot2::scale_fill_manual(values = colors_palette)
      # print(p)
    } else if (type == "hist") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(value, fill = variable)) +
        ggplot2::geom_histogram(alpha = 0.6, binwidth = 1, color = "black", show.legend = show_legend) + 
        ggplot2::theme_classic() + ggplot2::xlab("Coverage") +
        ggplot2::labs(fill = "Groups") +
        ggplot2::scale_fill_manual(values = colors_palette)
      # print(p)
    }
  }
  
  gc(verbose = FALSE)
  
  p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12,
                                                                       colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),
            axis.title.y = element_blank(), legend.title = element_blank())
  
  return(p + scMethrix_theme())
}

#--------------------------------------------------------------------------------------------------------------------------
#' Sparsity of sample
#'
#' @inheritParams generic_scMethrix_function
#' @param type Choose between 'box' (boxplot) or 'scatter' (scatterplot).
#' @return ggplot2 object
#' @examples
#' data('scMethrix_data')
#' plot_sparsity(scm = scMethrix_data)
#' @export
plot_sparsity <- function(scm = NULL, type = c("box", "scatter")) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  type <- match.arg(arg = type, choices = c("box", "scatter"), several.ok = FALSE)
  
  sparsity <- DelayedMatrixStats::colCounts(score(scm),value=NA)
  sparsity <- data.frame(variable = rownames(colData(scm)), Sparsity = sparsity/nrow(scm))
  
  if (type == "box") {
    p <- ggplot(sparsity, aes(x="", y=Sparsity)) + geom_boxplot() 
  }
  
  if (type == "scatter") {
    p <- ggplot(sparsity, aes(x="", y=Sparsity, color = variable)) + geom_point() 
  }

  p <- p + scMethrix_theme() + theme(axis.title.x = element_blank())

  return(p)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Plot descriptive statistics
#' @details plot descriptive statistics results from \code{\link{get_stats}}
#' @param plot_dat results from \code{\link{get_stats}}
#' @param what Can be \code{M} or \code{C}. Default \code{M}
#' @param stat Can be \code{mean} or \code{median}. Default \code{mean}
#' @param ignore_chr Chromsomes to ignore. Default \code{NULL}
#' @param samples Use only these samples. Default \code{NULL}
#' @param n_col number of columns. Passed to `facet_wrap`
#' @param n_row number of rows. Passed to `facet_wrap`
#' @return ggplot2 object
#' @seealso \code{\link{get_stats}}
#' @examples
#' data('scMethrix_data')
#' gs = get_stats(scMethrix_data)
#' plot_stats(gs)
#' @export
#'
plot_stats <- function(plot_dat, what = "Score", stat = "mean", ignore_chr = NULL,
                       samples = NULL, n_col = NULL, n_row = NULL) {
  
  if (is(plot_dat, "scMethrix")) plot_dat = get_stats(plot_dat)
  
  Chromosome <- . <- Sample_Name <- mean_meth <- sd_meth <- median_meth <- mean_cov <- sd_cov <- NULL
  median_cov <- measurement <- sd_low <- sd_high <- NULL
  stat <- match.arg(arg = stat, choices = c("mean", "median"))
  
  if ("Chr" %in% colnames(plot_dat)) {
    if (stat == "mean") {
      plot_dat[, which(grepl("^median", colnames(plot_dat))):=NULL]
    } else {
      plot_dat[, which(grepl("^mean", colnames(plot_dat))):=NULL]
    }
    
    if (!is.null(ignore_chr)) {
      plot_dat <- plot_dat[!Chromosome %in% ignore_chr]
    }
    
    if (!is.null(samples)) {
      plot_dat <- plot_dat[Sample_Name %in% samples]
    }
    
    colnames(plot_dat) <- c("Chromosome", "Sample_Name", "measurement",
                            "sd")
    plot_dat[, `:=`(measurement, as.numeric(as.character(measurement)))]
    plot_dat[, `:=`(sd, as.numeric(as.character(sd)))]
    plot_dat[, `:=`(sd_low, measurement - sd)]
    plot_dat[, `:=`(sd_high, measurement + sd)]
    plot_dat$sd_low <- ifelse(test = plot_dat$sd_low < 0, yes = 0,
                              no = plot_dat$sd_low)
    
    plot_dat_gg <- ggplot(data = plot_dat, aes(x = Chromosome, y = measurement)) +
      geom_errorbar(aes(ymin = sd_low, ymax = sd_high), col = "gray70") +
      geom_point(col = "maroon") + 
      facet_wrap(~Sample_Name, nrow = n_row, ncol = n_col) + 
      theme_minimal(base_size = 12) + 
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_text(hjust = 1, size = 10, colour = "black"),
            axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_blank())
  } else {
    if (stat == "mean") {
      plot_dat[, which(grepl("^median", colnames(plot_dat))):=NULL]
      plot_title <- paste("Mean",what)
    } else {
      plot_dat[, which(grepl("^mean", colnames(plot_dat))):=NULL]
      plot_title <- paste("Median",what)
    }
    
    colnames(plot_dat) <- c("Sample_Name", "measurement", "sd")
    plot_dat[, `:=`(measurement, as.numeric(as.character(measurement)))]
    plot_dat[, `:=`(sd, as.numeric(as.character(sd)))]
    plot_dat[, `:=`(sd_low, measurement - sd)]
    plot_dat[, `:=`(sd_high, measurement + sd)]
    plot_dat$sd_low <- ifelse(test = plot_dat$sd_low < 0, yes = 0,
                              no = plot_dat$sd_low)
    
    plot_dat_gg <- ggplot(data = plot_dat, aes(x = Sample_Name, y = measurement)) +
      geom_point(col = "maroon", size = 2) + 
      geom_errorbar(aes(ymin = sd_low, ymax = sd_high), col = "gray70") + 
      geom_point(col = "maroon") + theme_minimal(base_size = 12) + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"), 
            axis.text.y = element_text(size = 12, colour = "black"), 
            axis.title.y = element_blank()) + 
      ggtitle(label = plot_title)
  }
  
  return(plot_dat_gg  + scMethrix_theme())
}

plot_melissa <- function() {
  
  
  
}


plot_imap <- function(scm) {
  # 
  # x <- y <- NULL
  # 
  # umap <- get_matrix(scm,assay="umap")
  # 
  # df <- data.frame(x = scm$layout[,1],
  #                  y = scm$layout[,2])
  # 
  # ggplot(df, aes(x, y)) +
  #   geom_point()
  # 
}


#--------------------------------------------------------------------------------------------------------------------------
#' Plot dimensionality reduction
#' @param dim_red dimensionality reduction from an scMethrix object. Should be a matrix of two columns representing
#' the X and Y coordinates of the dim. red., with each row being a seperate sample
#' @param axis_labels A list of 'X' and 'Y' strings for labels, or NULL if no labels are desired
#' @param col_anno Column name of colData(m). Default NULL. Will be used as a factor to color different groups. Required \code{methrix} object
#' @param shape_anno Column name of colData(m). Default NULL. Will be used as a factor to shape different groups. Required \code{methrix} object
#' @param show_dp_labels Default FLASE
#' @return ggplot2 object
#' @importFrom graphics par mtext lines axis legend title
#' @examples
#' data('scMethrix_data')
#' scmpc = pca_scMethrix(scMethrix_data)
#' plot_pca(scmpc)
#' @export
plot_dim_red <- function(dim_red, col_anno = NULL, shape_anno = NULL, axis_labels = NULL, show_dp_labels = FALSE) {
  
  X <- Y <- color_me <- shape_me <- row_names <- NULL
  
  if (ncol(dim_red) != 2) {
    warning("More than two columns in the dimentionality reduction. Only the first two will be used")
    dim_red <- dim_red[,1:2]
  }
  
  if (!is.null(col_anno)) {
    col_anno_idx <- which(colnames(dim_red) == col_anno)
    if (length(col_anno_idx) == 0) {
      stop(paste0(col_anno, " not found in provided methrix object"))
    } else {
      colnames(dim_red)[col_anno_idx] <- "color_me"
    }
  }
  
  if (!is.null(shape_anno)) {
    shape_anno_idx <- which(colnames(dim_red) == shape_anno)
    if (length(shape_anno_idx) == 0) {
      stop(paste0(shape_anno, " not found in provided methrix object"))
    } else {
      colnames(dim_red)[shape_anno_idx] <- "shape_me"
    }
  }
  
  if (is.null(axis_labels)) {
    axis_labels = list(X="",Y="")
  }
  
  dim_red = as.data.frame(dim_red)
  colnames(dim_red) <- c("X", "Y")
  dim_red$row_names = rownames(dim_red)
  
  if (all(c("color_me", "shape_me") %in% colnames(dim_red))) {
    dimred_gg <- ggplot(data = dim_red, aes(x = X, y = Y, color = color_me,
                                            shape = shape_me, label = row_names)) + geom_point(size = 3) +
      labs(color = col_anno, shape = shape_anno) + scale_color_brewer(palette = "Dark2")
  } else if ("color_me" %in% colnames(dim_red)) {
    dimred_gg <- ggplot(data = dim_red, aes(x = X, y = Y, color = color_me,
                                            label = row_names)) + geom_point(size = 3) +
      labs(color = col_anno) + scale_color_brewer(palette = "Dark2")
  } else if ("shape_me" %in% colnames(dim_red)) {
    dimred_gg <- ggplot(data = dim_red, aes(x = X, y = Y, shape = shape_me,
                                            label = row_names)) + geom_point(size = 3) +
      labs(shape = shape_anno)
  } else {
    dimred_gg <- ggplot(data = as.data.frame(dim_red), aes(x = X, y = Y,
                                                           label = row_names)) + geom_point(size = 3, fill = "black",
                                                                                            color = "gray70")
  }
  
  dimred_gg <- dimred_gg  + theme_classic(base_size = 12) + xlab(axis_labels$X) + ylab(axis_labels$Y) + 
    theme(axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 12))
  
  if (show_dp_labels) {
    dimred_gg <- dimred_gg + geom_label(size = 4) 
  }
  
  return(dimred_gg)
  
}

#--------------------------------------------------------------------------------------------------------------------------
#' Plot PCA results
#' @inheritParams generic_scMethrix_function
#' @inheritParams plot_dim_red 
#' @param plot_vars Plot the variance explanation too
#' @param show_labels Show cell names on each data point. Default FLASE
#' @return ggplot2 object
#' @seealso [pca_scMethrix()] for dimensionality reduction
#' @importFrom graphics par mtext lines axis legend title barplot points
#' @examples
#' data('scMethrix_data')
#' scmpc = pca_scMethrix(scMethrix_data)
#' plot_pca(scmpc)
#' @export
plot_pca <- function(scm = NULL, col_anno = NULL, shape_anno = NULL, show_labels = FALSE, plot_vars = FALSE) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!("PCA" %in% reducedDimNames(scm))){
    stop("PCA results not present in scMethrix object. Run pca_scMethrix() first.")
  }
  
  dim_red = reducedDim(scm,type="PCA")[,1:2]
  pc_vars = scm@metadata$PCA_vars
  
  pc_x = colnames(dim_red)[1]
  pc_y = colnames(dim_red)[2]
  
  axis_labels = list(
    X = paste0(pc_x, " [", pc_vars[pc_x]*100, " %]"),
    Y = paste0(pc_y, " [", pc_vars[pc_y]*100, " %]")
  )
  
  pca_gg <- plot_dim_red(dim_red = dim_red, col_anno = col_anno, shape_anno = shape_anno, 
               show_dp_labels = show_labels, axis_labels = axis_labels)
  
  if (plot_vars) {
    par(mar = c(3, 4, 1, 4))
    b = barplot(height = pc_vars, names.arg = NA, col = "#2c3e50", ylim = c(0, 1), las = 2, axes = FALSE, ylab = "Variance Explained")
    points(x = b, y = cumsum(pc_vars), type = 'l', lty = 2, lwd = 1.2, xpd = TRUE, col = "#c0392b")
    points(x = b, y = cumsum(pc_vars), type = 'p', pch = 19, xpd = TRUE, col = "#c0392b")
    mtext(text = paste0("PC", 1:length(pc_vars)), side = 1, at = b, las = 2, line = 0.5, cex = 0.8)
    axis(side = 2, at = seq(0, 1, 0.1), line = 0, las = 2, cex.axis = 0.8)
    axis(side = 4, at = seq(0, 1, 0.1), line = 0, las = 2, cex.axis = 0.8)
    legend(x = "topleft", legend = "Cumulative", col = "#c0392b", pch = 19, lwd = 1, cex = 0.75, bty = "n")
  }
  
  return(pca_gg)
}

#' Plot tSNE results
#' @inheritParams plot_pca
#' @seealso [tsne_scMethrix()] for dimensionality reduction
#' @return ggplot2 object
#' @examples
#' data('scMethrix_data')
#' @export
plot_tsne <- function(scm = NULL, col_anno = NULL, shape_anno = NULL, show_labels = FALSE) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!("tSNE" %in% reducedDimNames(scm))){
    stop("t-SNE results not present in scMethrix object. Run tsne_scMethrix() first.")
  }
  
  plot_dim_red(dim_red = reducedDim(scm,type="tSNE"), col_anno = col_anno, shape_anno = shape_anno, 
               show_dp_labels = show_labels, axis_labels = NULL)
  
}

#' Plot UMAP results
#' @inheritParams plot_pca
#' @return ggplot2 object
#' @seealso [umap_scMethrix()] for dimensionality reduction
#' @examples
#' data('scMethrix_data')
#' @export
plot_umap <- function(scm = NULL, col_anno = NULL, shape_anno = NULL, show_labels = FALSE) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!("UMAP" %in% reducedDimNames(scm))){
    stop("UMAP results not present in scMethrix object. Run umap_scMethrix() first.")
  }

  plot_dim_red(dim_red = reducedDim(scm,type="UMAP"), col_anno = col_anno, shape_anno = shape_anno, 
               show_dp_labels = show_labels, axis_labels = NULL)
  
}

#------------------------------------------------------------------------------------------------------------
#' Evaluates imputations methods by NRMSE or AUC
#' @details Does stuff
#' @param sparse_prop numeric; A sparsity proportion between 0 and 1. E.g. 0.1 replaces 10% of the matrix with NA
#' @param imp_methods closure; The imputation methods to compare.
#' @param iterations integer; Number of iterations to test
#' @param type character; descriptive statistic. Can be either "AUC" or "RMSE". Default "RMSE"
#' @inheritParams generic_scMethrix_function
#' @return ggplot; The graph showing the NRMSE for each imputation method at each sparsity
#' @examples
#' data('scMethrix_data')
#' scMethrix_data <- impute_by_RF(scMethrix_data, new_assay="impute")
#' benchmark_imputation(scMethrix_data, assay="impute", sparse_prop = c(0.1,0.5,0.85))
#' @export
#' @import Metrics
benchmark_imputation <- function(scm = NULL, assay = "score", sparse_prop = seq(0.1, 0.9, 0.1), iterations = 3,
                                 imp_methods = c(iPCA = impute_by_iPCA, RF = impute_by_RF, kNN = impute_by_kNN),
                                 type = "RMSE") {
  
  results <- Sparsity <- NRMSE <- Imputation <- AUC <- NULL
  
  if (type == "AUC") {
    eq = Metrics::auc
  } else if (type == "RMSE") {
    eq = Metrics::rmse
  }
  
  
  for (prop in sparse_prop) {
    assays(scm)[["sparse"]] <- missForest::prodNA(get_matrix(scm,assay=assay),noNA=prop)
      
    for (i in 1:length(imp_methods)) {
      for (n in 1:iterations) {
          scm <- do.call(imp_methods[[i]], list(scm = scm,assay="sparse",new_assay="temp"))
          val <- do.call(eq,list(as.vector(get_matrix(scm,assay="temp")),
                                  as.vector(get_matrix(scm,assay=assay))))
          results = rbind(results,data.frame(Imputation=names(imp_methods)[i],Sparsity=prop,val=val))
    }}
  }

  results <- data.table(results)

  results = results[, .("mean" = mean(val), 
                         "sd" = sd(val)), 
                     by = c("Imputation", "Sparsity")]
  
  
  ggplot(results,aes(x=Sparsity, y=mean, color=Imputation)) +
    geom_line() + geom_point() + 
    geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) +
    xlab("Sparsity (proportion)") + ylab(type) + labs(fill = "Imputation") +
    scale_x_continuous(breaks=seq(0,1,.05)) + 
    scMethrix_theme() 
}

#------------------------------------------------------------------------------------------------------------
#' Theme for ggplot
#' @param base_size Size of text
#' @param base_family Family of text
#' @return ggplot element; data for the ggplot theme
#' @export
scMethrix_theme <- function(base_size = 12, base_family = "") {

  update_geom_defaults("line", list(size = 1.2))
  update_geom_defaults("point", list(size = 3))

  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      line =              element_line(colour = '#DADADA', size = 0.75, 
                                       linetype = 1, lineend = "butt"),
      rect =              element_rect(fill = "#F0F0F0", colour = "#F0F0F0", 
                                       size = 0.5, linetype = 1),
      text =              element_text(family = base_family, face = "plain",
                                       colour = "#656565", size = base_size,
                                       hjust = 0.5, vjust = 0.5, angle = 0, 
                                       lineheight = 0.9),
      plot.title =        element_text(size = rel(1.5), family = '' , 
                                       face = 'bold', hjust = -0.05, 
                                       vjust = 1.5, colour = '#3B3B3B'),
      axis.text =         element_text(),
      axis.ticks =        element_line(),
      axis.ticks.length.x = unit(.25, "cm"),
      axis.line =         element_line(colour = '#969696', size = 1, #TODO: Prefer only bottom line
                                       linetype = 1, lineend = "butt"),
      panel.grid.major =  element_line(colour = '#DADADA', size = 0.75, 
                                       linetype = 1, lineend = "butt"),
      panel.grid.minor =  element_blank(),
      plot.background =   element_rect(),
      panel.background =  element_rect(),
      legend.key =        element_rect(colour = '#DADADA'),
      complete = TRUE
    )
}