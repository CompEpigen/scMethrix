#' Format \code{\link{scMethrix}} matrix to long form data for plotting
#' 
#' @inheritParams generic_scMethrix_function
#' @param n_cpgs integer; Use these many random CpGs for plotting. Default 25000. Set it to \code{NULL} to use all - which can be memory expensive. The seed will be set to \code{n_cpgs} for consistency.
#' @param regions Granges; genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param pheno string; Col name of colData(m). Will be used as a factor to color different groups
#' @param na.rm boolean; remove NA values from the output
#' @return 'Long' matrix for methylation
#' @export
prepare_plot_data <- function(scm = NULL, assay="score", regions = NULL, n_cpgs = 25000, pheno = NULL, verbose = TRUE, na.rm = T){
  
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  .validateAssay(scm,assay)
  .validateType(regions,c("granges","null"))
  .validateType(n_cpgs,"integer")
  .validateType(pheno,c("string","null"))
  .validateType(na.rm,"boolean")

  Pheno <- Sample <- Value <- NULL
  
  #- Function code -----------------------------------------------------------------------------
  if (!is.null(regions)) {
    meth_sub <- subset_scMethrix(scm = scm, regions = regions)
    if (!is.null(n_cpgs)) {
      if(verbose) message("Randomly selecting ", n_cpgs, " sites")
      ids <- sample(x = seq_along(meth_sub), replace = FALSE, size = min(n_cpgs,
                                                                         nrow(meth_sub)))
      meth_sub <- get_matrix(scm = meth_sub[ids, ], assay = assay, add_loci = FALSE)
    } else {
      meth_sub <- get_matrix(scm = meth_sub, assay = assay, add_loci = FALSE)
    }
  } else if (!is.null(n_cpgs)) {
    if(verbose) message("Randomly selecting ", n_cpgs, " sites")
    
    set.seed(n_cpgs)
    ids <- sample(x = seq_along(scm), replace = FALSE, size = min(n_cpgs,
                                                                nrow(scm)))
    meth_sub <- get_matrix(scm = scm[ids, ], assay = assay, add_loci = FALSE)
  } else {
    meth_sub <- get_matrix(scm = scm, assay = assay, add_loci = FALSE)
  }
  
  meth_sub <- as.data.frame(meth_sub)
  data.table::setDT(x = meth_sub)
  plot.data <- suppressWarnings(data.table::melt(meth_sub))
  
  colnames(plot.data) <- c("Sample", "Value")
  
  if (!is.null(pheno)) {
    if (pheno %in% colnames(colData(scm))) {
      #colnames(meth_sub) <- as.character(scm@colData[, pheno])
      #TODO: make sure the order is correct
      plot.data[, `:=`(Pheno, rep(colData(scm)[,pheno],each=nrow(meth_sub)))]
      plot.data[, `:=`(Pheno, as.factor(Pheno))]
    } else {
      stop("Please provide a valid phenotype annotation column.")
    }
  } else {
    plot.data[, `:=`(Pheno, Sample)]
  }

  if(na.rm) plot.data <- plot.data[!is.na(Value)]
  
  gc(verbose = FALSE)
  return(plot.data)
  
}

#' Getter for plot palette colors
#' @param n_row Number of colors
#' @param col_palette String for RColorBrewer palette name  
#' @return RColorBrewer palette
get_palette <- function(n_row, col_palette){
  
  #- Input Validation --------------------------------------------------------------------------
  .validateType(n_row,"integer")
  .validateType(col_palette,"string")
  
  if (n_row == 0) {
    stop("Zero colors present in the palette")
  }
  
  if (!col_palette %in% row.names(RColorBrewer::brewer.pal.info)){
    stop("Please provide a valid RColorBrewer palettte. Possible values are: ", paste0(row.names(RColorBrewer::brewer.pal.info)), sep=", ")
  }
  
  #- Function code -----------------------------------------------------------------------------
  color_pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[col_palette,
                                                                                                  "maxcolors"], col_palette))(n_row)
  return(color_pal)
}

#' Getter for plot shapes. Shapes selected for optimal distinction and taken from:
#' @details http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
#' @param n_row Number of shapes. Max of 15.
#' @return list of shapes (by integer)
get_shape <- function(n_row) {
  .validateType(n_row,"integer")
  shapes <- c(21:25,3,4,7:14)
  return(shapes[1:n_row])
}



#' Violin Plot for \eqn{\beta}-Values
#' @inheritParams prepare_plot_data
#' @param col_palette string; Name of the RColorBrewer palette to use for plotting.
#' @param show_legend boolean; Display the legend on the plot
#' @return ggplot2 object
#' @export
#' @import ggplot2
#' @examples
#' data('scMethrix_data')
#' plot_violin(scm = scMethrix_data)
plot_violin <- function(scm = NULL, assay="score", regions = NULL, n_cpgs = 25000, pheno = NULL,
                        col_palette = "RdYlGn", show_legend = TRUE, verbose = TRUE) {
  
  #- Input Validation --------------------------------------------------------------------------
  Sample <- Value <- Pheno <- NULL
  
  .validateExp(scm)
  .validateAssay(scm,assay)
  .validateType(regions,c("granges","null"))
  .validateType(n_cpgs,"integer")
  .validateType(pheno,c("string","null"))
  .validateType(col_palette,"string")
  .validateType(show_legend,"boolean")

  if (is.null(regions)) regions = rowRanges(scm)
  
  #- Function code -----------------------------------------------------------------------------
  plot.data <- prepare_plot_data(scm=scm, assay = assay, regions = regions, n_cpgs = n_cpgs, pheno = pheno)
  
  col_palette <- get_palette(ncol(scm), col_palette)
  # generate the violin plot
  
  p <- ggplot2::ggplot(plot.data, ggplot2::aes(x = Sample, y = Value, fill = Pheno)) + 
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
plot_density <- function(scm = NULL, assay = "score", regions = NULL, n_cpgs = 25000, pheno = NULL,
                         col_palette = "RdYlGn", show_legend = TRUE, verbose = TRUE) {
  
  #- Input Validation --------------------------------------------------------------------------
  Value <- Pheno <- NULL
  
  .validateExp(scm)
  .validateAssay(scm,assay)
  .validateType(regions,c("granges","null"))
  .validateType(n_cpgs,"integer")
  .validateType(pheno,c("string","null"))
  .validateType(col_palette,"string")
  .validateType(show_legend,"boolean")
  
  if (is.null(regions)) regions = rowRanges(scm)
  
  #- Function code -----------------------------------------------------------------------------
  plot.data <- prepare_plot_data(scm=scm, regions = regions, n_cpgs = n_cpgs, pheno = pheno)
  col_palette <- get_palette(ncol(scm), col_palette)
  
  # generate the density plot

  p <- ggplot2::ggplot(plot.data, ggplot2::aes(Value, color = Pheno)) + geom_density(lwd = 1, position = "identity", show.legend = show_legend) + ggplot2::theme_classic() +
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
#' @param max_cov integer; Maximum coverage value to be plotted.
#' @param type string; Choose between 'histogram' (histogram) or 'density' (density plot).
#' @param obs_lim integer; The maximum number of observations (sites*samples) to use. If the dataset is larger that this,
#' random sites will be selected from the genome.
#' @return ggplot2 object
#' @examples
#' data('scMethrix_data')
#' plot_coverage(scm = scMethrix_data)
#' @export
plot_coverage <- function(scm = NULL, type = c("histogram", "density"), pheno = NULL,
                          max_cov = 100, obs_lim = 1e+06, col_palette = "RdYlGn", show_legend = TRUE, verbose = TRUE) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  type <- .validateArg(type, plot_coverage)
  .validateType(pheno,c("string","null"))
  .validateType(max_cov,"integer")
  .validateType(obs_lim,"integer")
  .validateType(col_palette,"string")
  .validateType(show_legend,"boolean")
  
  Value <- Sample <- Pheno <- NULL
  
  colors_palette <- get_palette(ncol(scm), col_palette)
  
  #- Function code -----------------------------------------------------------------------------
  if (matrixStats::product(dim(scm)) > obs_lim) {
    message("The dataset is bigger than the size limit. A random subset of the object will be used that contains ~",
            obs_lim, " observations.")
    n_rows <- trunc(obs_lim/ncol(scm))
    sel_rows <- sample(seq_len(nrow(scm)), size = n_rows,
                       replace = FALSE)
  } else {
    sel_rows <- seq_len(nrow(scm))
  }

  plot.data <- prepare_plot_data(scm = scm[sel_rows, ], assay = "counts", pheno = pheno, na.rm = F)
  setnafill(plot.data,fill=0,cols="Value")
  
  plot.data <- plot.data[Value <= max_cov, ]

  # generate the plots
  if (is.null(pheno)) {
    if (type == "density") {
      p <- ggplot2::ggplot(plot.data, aes(Value, color = Sample)) +
        ggplot2::geom_density(alpha = 0.5, adjust = 1.5, lwd = 1, show.legend = show_legend,
                              position = "identity") + ggplot2::theme_classic() + ggplot2::xlab("Coverage") +
        ggplot2::scale_fill_manual(values = colors_palette)
      
    } else if (type == "histogram") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(Value, fill = Sample)) + 
        ggplot2::geom_histogram(alpha = 0.6, binwidth = 1, color = "black", show.legend = show_legend) + 
        ggplot2::theme_classic() +
        ggplot2::xlab("Coverage")+
        ggplot2::scale_fill_manual(values = colors_palette)
      # print(p)
    }
  } else {
    if (type == "density") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(Value, color = Pheno)) +
        ggplot2::geom_density(alpha = 0.6, adjust = 1.5, lwd = 1, show.legend = show_legend,
                              position = "identity") + ggplot2::theme_classic() + ggplot2::xlab("Coverage") +
        ggplot2::scale_fill_manual(values = colors_palette)
      # print(p)
    } else if (type == "histogram") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(Value, fill = Pheno)) +
        ggplot2::geom_histogram(alpha = 0.6, binwidth = 1, color = "black", show.legend = show_legend) + 
        ggplot2::theme_classic() + ggplot2::xlab("Coverage") +
        ggplot2::scale_fill_manual(values = colors_palette)
      # print(p)
    }
  }
  
  gc(verbose = FALSE)
  
  p <- p + ggplot2::theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12,
                                                                       colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),
            axis.title.y = element_blank(), legend.title = element_blank())
  
  return(p + scMethrix_theme())
}

#--------------------------------------------------------------------------------------------------------------------------
#' Sparsity of sample
#' inheritParams generic_plot_function
#' @inheritParams plot_violin
#' @param type string; Choose between 'box' (boxplot) or 'scatter' (scatterplot).
#' @return ggplot2 object
#' @examples
#' data('scMethrix_data')
#' plot_sparsity(scm = scMethrix_data)
#' @export
plot_sparsity <- function(scm = NULL, assay = "score", type = c("box", "scatter"), pheno = NULL, verbose = TRUE) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  type <- .validateArg(type,plot_sparsity)
  .validateAssay(scm,assay)
  .validateType(pheno,c("string","null"))
  
  Sparsity <- variable <- NULL
  
  sparsity <- DelayedMatrixStats::colCounts(get_matrix(scm,assay=assay),value=NA)
  
  #- Function code -----------------------------------------------------------------------------
  if (!is.null(pheno)) {
    if (pheno %in% colnames(colData(scm))) {
      pheno <- as.character(scm@colData[, pheno])
      sparsity <- data.frame(Phenotype = pheno, Sparsity = sparsity/nrow(scm))
      p <- ggplot2::ggplot(sparsity, aes(x=pheno, y=Sparsity, color = pheno))
    } else {
      stop("Please provide a valid phenotype annotation column.")
    }
  } else {
    sparsity <- data.frame(Sparsity = sparsity/nrow(scm))
    p <- ggplot2::ggplot(sparsity, aes(x="", y=Sparsity))
  }

  if (type == "box") {p <- p + ggplot2::geom_boxplot()
  } else if (type == "scatter") {p <- p + ggplot2::geom_point() }

  p <- p + scMethrix_theme() + ggplot2::theme(axis.title.x = element_blank())

  return(p)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Plot descriptive statistics
#' @details plot descriptive statistics results from \code{\link{get_stats}}
#' @inheritParams plot_violin
#' @param scm scMethrix; \code{\link{get_stats}} will be run for the specified assay
#' @param stat string; Can be \code{mean} or \code{median}. Default \code{mean}
#' @param ignore_chr string; Chromsomes to ignore. If NULL, all chromosome will be used. Default \code{NULL}
#' @param ignore_samples list of strings; Samples to ignore.  If NULL, all samples will be used. Default \code{NULL}
#' @param n_col integer; number of columns. Passed to `facet_wrap`
#' @param n_row integer; number of rows. Passed to `facet_wrap`
#' @param per_chr boolean; plot per chromosome
#' @return ggplot2 object
#' @seealso \code{\link{get_stats}}
#' @examples
#' data('scMethrix_data')
#' plot_stats(scMethrix_data)
#' @export
#'
plot_stats <- function(scm, assay = "score", stat = c("mean", "median"), per_chr = FALSE, ignore_chr = NULL,
                       ignore_samples = NULL, n_col = NULL, n_row = NULL, pheno = NULL, verbose = TRUE) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  stat <- .validateArg(stat,plot_stats)
  .validateType(per_chr,"boolean")
  .validateType(ignore_chr,c("string","null"))
  .validateType(ignore_samples,c("string","null"))
  .validateType(n_col,c("integer","null"))
  .validateType(n_row,c("integer","null"))
  
  Chromosome <- . <- Sample <- mean_meth <- sd_meth <- median_meth <- mean_cov <- sd_cov <- NULL
  median_cov <- measurement <- sd_low <- sd_high <- NULL

  plot_dat = get_stats(scm, assay = assay, per_chr = per_chr, ignore_chr = ignore_chr, ignore_samples = ignore_samples)

  #- Function code -----------------------------------------------------------------------------
  if (per_chr) {
    if (stat == "mean") {
      plot_dat[, which(grepl("^median", colnames(plot_dat))):=NULL]
    } else {
      plot_dat[, which(grepl("^mean", colnames(plot_dat))):=NULL]
    }

    colnames(plot_dat) <- c("Chromosome", "Sample", "measurement",
                            "sd")
    plot_dat[, `:=`(measurement, as.numeric(as.character(measurement)))]
    plot_dat[, `:=`(sd, as.numeric(as.character(sd)))]
    plot_dat[, `:=`(sd_low, measurement - sd)]
    plot_dat[, `:=`(sd_high, measurement + sd)]
    plot_dat$sd_low <- ifelse(test = plot_dat$sd_low < 0, yes = 0,
                              no = plot_dat$sd_low)
    
    plot_dat_gg <- ggplot(data = plot_dat, aes(x = Chromosome, y = measurement)) +
      ggplot2::geom_errorbar(aes(ymin = sd_low, ymax = sd_high), col = "gray70") +
      ggplot2::geom_point(col = "maroon") + 
      ggplot2::facet_wrap(~Sample, nrow = n_row, ncol = n_col) + 
      ggplot2::theme_minimal(base_size = 12) + 
      ggplot2::theme(axis.title.x = element_blank(), 
            axis.text.x = element_text(hjust = 1, size = 10, colour = "black"),
            axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_blank())
  } else {
    if (stat == "mean") {
      plot_dat[, which(grepl("^median", colnames(plot_dat))):=NULL]
      plot_title <- paste("Mean",assay)
    } else {
      plot_dat[, which(grepl("^mean", colnames(plot_dat))):=NULL]
      plot_title <- paste("Median",assay)
    }
    
    colnames(plot_dat) <- c("Sample", "measurement", "sd")
    plot_dat[, `:=`(measurement, as.numeric(as.character(measurement)))]
    plot_dat[, `:=`(sd, as.numeric(as.character(sd)))]
    plot_dat[, `:=`(sd_low, measurement - sd)]
    plot_dat[, `:=`(sd_high, measurement + sd)]
    plot_dat$sd_low <- ifelse(test = plot_dat$sd_low < 0, yes = 0,
                              no = plot_dat$sd_low)
    
    plot_dat_gg <- ggplot2::ggplot(data = plot_dat, aes(x = Sample, y = measurement)) +
      ggplot2::geom_point(col = "maroon", size = 2) + 
      ggplot2::geom_errorbar(aes(ymin = sd_low, ymax = sd_high), col = "gray70") + 
      ggplot2::geom_point(col = "maroon") + theme_minimal(base_size = 12) + 
      ggplot2::theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"), 
            axis.text.y = element_text(size = 12, colour = "black"), 
            axis.title.y = element_blank()) + 
      ggplot2::ggtitle(label = plot_title)
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
#' @inheritParams generic_scMethrix_function
#' @param dim_red string; name of adimensionality reduction from an scMethrix object. Should be a matrix of two columns representing
#' the X and Y coordinates of the dim. red., with each row being a seperate sample
#' @param axis_labels list of strings; A list of 'X' and 'Y' strings for labels, or NULL if no labels are desired
#' @param color_anno string; Column name of colData(m). Default NULL. Will be used as a factor to color different groups. Required \code{methrix} object
#' @param shape_anno string; Column name of colData(m). Default NULL. Will be used as a factor to shape different groups. Required \code{methrix} object
#' @param show_dp_labels boolean; Flag to show the labels for dots. Default FALSE
#' @return ggplot2 object
#' @importFrom graphics par mtext lines axis legend title
#' @export
plot_dim_red <- function(scm, dim_red, color_anno = NULL, shape_anno = NULL, axis_labels = NULL, show_dp_labels = FALSE, verbose = TRUE) {

  #- Input Validation --------------------------------------------------------------------------
  X <- Y <- Color <- Shape <- color <- shape <- Sample <- row_names <- NULL
  
  .validateExp(scm)
  .validateType(dim_red,"string")
  if (!(dim_red %in% reducedDimNames(scm))) stop("Invalid dim_red specified. '",dim_red,"' does not exist in the experiment.")
  .validateType(color_anno,c("string","null"))
  .validateType(shape_anno,c("string","null"))
  .validateType(unlist(axis_labels),c("string","null"))
  .validateType(show_dp_labels,"boolean")
  
  dim_red <- reducedDim(scm,type=dim_red)
  
  if (ncol(dim_red) != 2) {
    warning("More than two columns in the dimentionality reduction. Only the first two will be used")
    dim_red <- dim_red[,1:2]
  }
  
  dim_red = as.data.frame(dim_red)
  colnames(dim_red) <- c("X", "Y")
  dim_red$Sample = rownames(dim_red)
  
  #- Function code -----------------------------------------------------------------------------
  if (!is.null(color_anno)) {
    if (color_anno  %in% colnames(colData(scm))) {
      dim_red$Color <- as.factor(unlist(as.data.table(colData(scm))[,color_anno, with=FALSE])) #TODO: make colData a data.table
      colors <- scale_color_manual(values= get_palette(length(unique(dim_red$Color)),col_palette = "Dark2"))
    } else {
      stop(paste0(color_anno, " not found in provided scMethrix object"))
    }
  }
    
  if (!is.null(shape_anno)) {
    if (shape_anno %in% colnames(colData(scm))) {
      dim_red$Shape <- as.factor(unlist(as.data.table(colData(scm))[,shape_anno, with=FALSE])) #TODO: make colData a data.table
      shapes <- scale_shape_manual(values = get_shape(length(unique(dim_red$Shape))))
    } else {
      stop(paste0(shape_anno, " not found in provided scMethrix object"))
    }
  }  

  if (is.null(axis_labels)) {
    axis_labels = list(X="",Y="")
  }
  
  if (all(c("Color", "Shape") %in% colnames(dim_red))) {
    dimred_gg <- ggplot2::ggplot(data = dim_red, aes(x = X, y = Y, color = Color,
                                            shape = Shape, label = Sample)) + geom_point(size = 3) +
      labs(color = color_anno, shape = shape_anno) + scale_color_brewer(palette = "Dark2")
  } else if ("Color" %in% colnames(dim_red)) {
    dimred_gg <- ggplot2::ggplot(data = dim_red, aes(x = X, y = Y, color = Color,
                                            label = Sample)) + geom_point(size = 3) +
      labs(color = color_anno) + scale_color_brewer(palette = "Dark2")
  } else if ("Shape" %in% colnames(dim_red)) {
    dimred_gg <- ggplot2::ggplot(data = dim_red, aes(x = X, y = Y, shape = Shape,
                                            label = Sample)) + geom_point(size = 3) +
      labs(shape = shape_anno)
  } else {
    dimred_gg <- ggplot2::ggplot(data = as.data.frame(dim_red), aes(x = X, y = Y,
                                                           label = Sample)) + geom_point(size = 3, fill = "black",
                                                                                            color = "gray70")
  }
  
  dimred_gg <- dimred_gg  + ggplot2::theme_classic(base_size = 12) + 
    ggplot2::xlab(axis_labels$X) + ggplot2::ylab(axis_labels$Y) + 
    ggplot2::theme(axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 12))
  
  if (show_dp_labels) {
    dimred_gg <- dimred_gg + ggplot2::geom_label(size = 4) 
  }
  
  dimred_gg <- dimred_gg + color + shape
  
  return(dimred_gg)
  
}

#--------------------------------------------------------------------------------------------------------------------------
# Plot PCA results
# @inheritParams generic_scMethrix_function
# @inheritParams plot_dim_red 
# @param plot_vars Plot the variance explanation too
# @param show_labels Show cell names on each data point. Default FLASE
# @return ggplot2 object
# @seealso [pca_scMethrix()] for dimensionality reduction
# @importFrom graphics par mtext lines axis legend title barplot points
# @examples
# data('scMethrix_data')
# scmpc = dim_red_scMethrix(scMethrix_data,type="PCA")
# plot_pca(scmpc)
# @export
# plot_pca <- function(scm = NULL, col_anno = NULL, shape_anno = NULL, show_labels = FALSE, plot_vars = FALSE) {
#   
#   if (!is(scm, "scMethrix")){
#     stop("A valid scMethrix object needs to be supplied.")
#   }
#   
#   if (!("PCA" %in% reducedDimNames(scm))){
#     stop("PCA results not present in scMethrix object. Run pca_scMethrix() first.")
#   }
#   
#   dim_red = reducedDim(scm,type="PCA")[,1:2]
#   pc_vars = scm@metadata$PCA_vars
#   
#   pc_x = colnames(dim_red)[1]
#   pc_y = colnames(dim_red)[2]
#   
#   axis_labels = list(
#     X = paste0(pc_x, " [", pc_vars[pc_x]*100, " %]"),
#     Y = paste0(pc_y, " [", pc_vars[pc_y]*100, " %]")
#   )
#   
#   pca_gg <- plot_dim_red(scm, dim_red = "PCA", col_anno = col_anno, shape_anno = shape_anno, 
#                show_dp_labels = show_labels, axis_labels = axis_labels)
#   
#   if (plot_vars) {
#     par(mar = c(3, 4, 1, 4))
#     b = barplot(height = pc_vars, names.arg = NA, col = "#2c3e50", ylim = c(0, 1), las = 2, axes = FALSE, ylab = "Variance Explained")
#     points(x = b, y = cumsum(pc_vars), type = 'l', lty = 2, lwd = 1.2, xpd = TRUE, col = "#c0392b")
#     points(x = b, y = cumsum(pc_vars), type = 'p', pch = 19, xpd = TRUE, col = "#c0392b")
#     mtext(text = paste0("PC", 1:length(pc_vars)), side = 1, at = b, las = 2, line = 0.5, cex = 0.8)
#     axis(side = 2, at = seq(0, 1, 0.1), line = 0, las = 2, cex.axis = 0.8)
#     axis(side = 4, at = seq(0, 1, 0.1), line = 0, las = 2, cex.axis = 0.8)
#     legend(x = "topleft", legend = "Cumulative", col = "#c0392b", pch = 19, lwd = 1, cex = 0.75, bty = "n")
#   }
#   
#   return(pca_gg)
# }

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
#' \dontrun{
#' scMethrix_data <- impute_regions(scMethrix_data, new_assay="impute",type="RF")
#' benchmark_imputation(scMethrix_data, assay="impute", sparse_prop = c(0.1,0.5,0.85))
#' }
#' @export
#' @import Metrics
benchmark_imputation <- function(scm = NULL, assay = "score", sparse_prop = seq(0.1, 0.9, 0.1), iterations = 3,
                                 imp_methods = c(iPCA = function(...) impute_regions(type="iPCA",...), 
                                                 RF = function(...) impute_regions(type="RF",...), 
                                                 kNN = function(...) impute_regions(type="kNN",...)),
                                 type = "RMSE") {
  
  #- Input Validation --------------------------------------------------------------------------
  . <- results <- Sparsity <- NRMSE <- Imputation <- AUC <- NULL
  
  #- Function code -----------------------------------------------------------------------------
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
  
  
  ggplot2::ggplot(results,aes(x=Sparsity, y=mean, color=Imputation)) +
    ggplot2::geom_line() + ggplot2::geom_point() + 
    ggplot2::geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) +
    ggplot2::xlab("Sparsity (proportion)") + ggplot2::ylab(type) + ggplot2::labs(fill = "Imputation") +
    ggplot2::scale_x_continuous(breaks=seq(0,1,.05)) + 
    scMethrix_theme() 
}

#------------------------------------------------------------------------------------------------------------
#' Theme for ggplot
#' @param base_size integer; Size of text
#' @param base_family string; Family of text
#' @return ggplot element; data for the ggplot theme
#' @export
scMethrix_theme <- function(base_size = 12, base_family = "") {

  update_geom_defaults("line", list(size = 1.2))
  update_geom_defaults("point", list(size = 3))

  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
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
