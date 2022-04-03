#---- prepare_plot_data ------------------------------------------------------------------------------------------------
#' Formats an [`scMethrix`] matrix to long form data for plotting
#' @inheritParams generic_scMethrix_function
#' @param n_cpgs `integer`; Use these many random CpGs for plotting. Default = `25000`. Set it to `NULL` to use all - which can be memory expensive. The seed will be set to `n_cpgs` for consistency.
#' @param pheno `string`; col name of `colData(scm)`. Will be used as a factor to color different groups
#' @param na.rm `boolean`; remove NA values from the output
#' @return 'Long' matrix for methylation
#' @export
prepare_plot_data <- function(scm = NULL, assay="score", n_cpgs = 25000, pheno = NULL, verbose = TRUE, na.rm = T){
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  .validateAssay(scm,assay)
  .validateType(n_cpgs,"integer")
  .validateType(pheno,c("string","null"))
  .validateType(na.rm,"boolean")

  Pheno <- Sample <- Value <- NULL
  
  if (n_cpgs > nrow(scm)) n_cpgs = NULL
  
  #---- Function code ------------------------------------------------------
  if (!is.null(n_cpgs)) {
    if(verbose) message("Randomly selecting ", n_cpgs, " sites")
    set.seed(n_cpgs)
    ids <- sort(sample(x = seq_along(scm), replace = FALSE, size = min(n_cpgs, nrow(scm))))
    meth_sub <- get_matrix(scm = scm[ids, ], assay = assay, add_loci = FALSE)
  } else {
    meth_sub <- get_matrix(scm = scm, assay = assay, add_loci = FALSE)
  }
  
  meth_sub <- as.data.table(meth_sub)
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

#---- .getPalette ------------------------------------------------------------------------------------------------------
#' Generates a palette of RGB colors via [`colorspace`](https://colorspace.r-forge.r-project.org//index.html).
#' @description 
#' This function wraps the color generation functions from [`colorspace`](https://cran.r-project.org/web/packages/colorspace/index.html): [`sequential_hcl()`][colorspace::sequential_hcl()], [`diverging_hcl()`][colorspace::diverging_hcl()], and [`qualitative_hcl()`][colorspace::qualitative_hcl()]. For simplicity, as the palette IDs are unique, this function will automatically select the appropriate function to generate the palettes.
#' @details
#' All [`colorspace`](https://colorspace.r-forge.r-project.org//index.html) palettes can be visualized with:
#' ```
#' library("colorspace")
#' hcl_palettes(plot = TRUE)
#' ```
#' 
#' Custom palettes can also be [added](https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html#registering-your-own-palettes-1).
#' 
#' To easily visualize the colors, use `scales::show_col(.getPalette())`
#' 
#' @param nColors `integer`; Number of colors. Default = `10`.
#' @param paletteID `string`; the ID of the [`colorspace`](https://colorspace.r-forge.r-project.org//index.html) palette. Default = `"Dark Mint"`.
#' @param ... Additional arguments for palette generation.
#' @return vector of HEX colors
#' @examples 
#' palette <- .getPalette(nColors = 3, paletteID = "Plasma")
#' scales::show_col(palette)#' 
#' df <- mtcars[, c("mpg", "cyl", "wt")]
#' df$cyl <- as.factor(df$cyl)
#' 
#' ggplot2::ggplot(df, aes(x=wt, y=mpg, group=cyl)) + geom_point()
#' ggplot2::ggplot(df, aes(x=wt, y=mpg, group=cyl)) + geom_point(aes(color=cyl)) + scale_color_manual(values= palette)
#' @noRd
.getPalette <- function(nColors = 10, paletteID = "Dark Mint", ...){
  
  #---- Input validation ---------------------------------------------------
  .validatePackageInstall("colorspace")
  
  .validateType(paletteID,"string")
  .validateType(nColors,c("integer"))
  .validateValue(nColors,"> 1")

  validPalettes <- row.names(colorspace::hcl_palettes())
  if (!(paletteID %in% validPalettes)) {
    stop("Invalid palette. Must be one of: '", paste(validPalettes,collapse = "', '"),"'.")
  }
  
  #---- Function code ------------------------------------------------------
  hcl <- colorspace::hcl_palettes(palette = paletteID)
  hcl <- as.character(hcl$type)
  hcl <- gsub( " .*$", "", hcl )
  
  if (hcl == "Sequential")        palette <- colorspace::sequential_hcl (nColors, palette = paletteID, ...) 
  else if (hcl == "Diverging")    palette <- colorspace::diverging_hcl  (nColors, palette = paletteID, ...) 
  else if (hcl == "Qualitative")  palette <- colorspace::qualitative_hcl(nColors, palette = paletteID, ...) 
  
  return(palette)
}

#---- .getShapes -------------------------------------------------------------------------------------------------------
#' Generates list of optimally distinct pre-selected shapes.
#' @description For general list of shapes, see [here](http://sape.inf.usi.ch/quick-reference/ggplot2/shape).
#' @details
#' The order of the shapes (left-right, top-bottom):
#' ```
#' shapes <- .getShapes()
#' d <- data.frame(p = shapes, i = 1:length(shapes)-1)
#' ggplot() +
#'  scale_y_continuous(name="") +
#'  scale_x_continuous(name="") +
#'  scale_y_reverse() +
#'  scale_shape_identity() +
#'  geom_point(data=d, mapping=aes(x=i%%16, y=i%/%16, shape=p), size=5, fill="red") +
#'  geom_text(data=d, mapping=aes(x=i%%16, y=i%/%16+0.25, label=p), size=3) +
#'  theme_void()
#'```
#' @param nShapes `integer`; Number of shapes. Max of 57.
#' @return `list(integer)`; the specified number of shapes. If nShapes = `NULL`, returns all possible shapes.
#' @examples 
#' shapes <- .getShapes(3) 
#' df <- mtcars[, c("mpg", "cyl", "wt")]
#' df$cyl <- as.factor(df$cyl)
#' 
#' ggplot2::ggplot(df, aes(x=wt, y=mpg, group=cyl)) + geom_point()
#' ggplot2::ggplot(df, aes(x=wt, y=mpg, group=cyl)) + geom_point(aes(shape=cyl)) + scale_shape_manual(values= shapes)
#' @noRd
.getShapes <- function(nShapes = NULL) {
  
  #---- Input validation ---------------------------------------------------
  .validateType(nShapes,c("integer","null"))
  
  #---- Function code ------------------------------------------------------
  shapes <- c(15:25,3,4,7:14,97:122,48:57)
  
  if (is.null(nShapes)) return (shapes)
  
  if (length(nShapes) > length(shapes)) stop("Invalid shapes. Max of ",length(shapes)," shapes possible", call. = FALSE)
  
  return(shapes[1:nShapes])
}

#---- plot_violin ------------------------------------------------------------------------------------------------------
#' Violin Plot for `\beta`-Values
#' @inheritParams prepare_plot_data
#' @param n_cpgs `integer`; The number of CpGs to for plotting. Default = `25000`.
#' @param paletteID `string`; Name of the `colorspace` palette to use for plotting.
#' @param show_legend `boolean`; Display the legend on the plot
#' @param ... Additional parameters to feed to [scMethrix_theme()]
#' @return [`ggplot2::ggplot2`] object
#' @export
#' @import ggplot2
#' @examples
#' data('scMethrix_data')
#' plot_violin(scm = scMethrix_data)
plot_violin <- function(scm = NULL, assay="score", n_cpgs = 25000, pheno = NULL,
                        paletteID = "Dark Mint", show_legend = FALSE, verbose = TRUE,...) {
  
  #---- Input validation ---------------------------------------------------
  Sample <- Value <- Pheno <- NULL
  
  .validateExp(scm)
  .validateAssay(scm,assay)
  .validateType(n_cpgs,"integer")
  .validateType(pheno,c("string","null"))
  .validateType(paletteID,"string")
  .validateType(show_legend,"boolean")

  #---- Function code ------------------------------------------------------
  plot.data <- prepare_plot_data(scm=scm, assay = assay, n_cpgs = n_cpgs, pheno = pheno)
  
  palette <- .getPalette(ncol(scm), paletteID)
  # generate the violin plot
  
  p <- ggplot2::ggplot(plot.data, ggplot2::aes(x = Sample, y = Value, fill = Pheno)) + 
    ggplot2::geom_violin(alpha = 0.8, show.legend = show_legend) + ggplot2::theme_classic(base_size = 14) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::xlab(pheno) + ggplot2::ylab(expression(beta * "-value")) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12,
                                                                     colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),
          axis.title.y = element_blank(), legend.title = element_blank())
  
  return(p + scMethrix_theme(...))
}

#---- plot_density -----------------------------------------------------------------------------------------------------
#' Density Plot of `\beta`-Values
#'
#' @inheritParams plot_violin
#' @param na.rm `boolean`; Remove NA values from the plot
#' @return [`ggplot2::ggplot2`] object
#' @export
#' @examples
#' data('scMethrix_data')
#' plot_density(scm = scMethrix_data)
plot_density <- function(scm = NULL, assay = "score", n_cpgs = 25000, pheno = NULL,
                         paletteID = "Dark Mint", show_legend = FALSE, verbose = TRUE, na.rm = T,...) {
  
  #---- Input validation ---------------------------------------------------
  Value <- Pheno <- NULL
  
  .validateExp(scm)
  .validateAssay(scm,assay)
  .validateType(n_cpgs,"integer")
  .validateType(pheno,c("string","null"))
  .validateType(paletteID,"string")
  .validateType(show_legend,"boolean")
  
  #---- Function code ------------------------------------------------------
  plot.data <- prepare_plot_data(scm=scm, assay = assay, n_cpgs = n_cpgs, pheno = pheno)
  palette <- .getPalette(ncol(scm), paletteID)

    # generate the density plot

  p <- ggplot2::ggplot(plot.data, ggplot2::aes(Value, color = Pheno)) + geom_density(lwd = 1, position = "identity", show.legend = show_legend,kernel="cosine",na.rm = na.rm) + ggplot2::theme_classic() +
    ggplot2::xlab("Methylation") + ggplot2::ylab("Density") + ggplot2::theme_classic(base_size = 14) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::xlab(expression(beta * "-value")) + theme(axis.title.x = element_blank(), 
                                     axis.text.x = element_text(size = 12, colour = "black"), 
                                     axis.text.y = element_text(size = 12, colour = "black"), 
                                     axis.title.y = element_blank(), legend.title = element_blank())+
                                     ggplot2::scale_color_manual(values = .getPalette(length(levels(plot.data$Pheno)),paletteID = paletteID))
    
  
  gc(verbose = FALSE)
  
  return(p + scMethrix_theme(...))
}

#---- plot_coverage ----------------------------------------------------------------------------------------------------
#' Coverage QC Plots
#' @inheritParams plot_violin
#' @param max_cov `integer`; Maximum coverage value to be plotted.
#' @param type `string`; Choose between `histogram` or `density` plot.
#' @param obs_lim `integer`; The maximum number of observations (`sites*samples`) to use. If the dataset is larger that this, random sites will be selected from the genome.
#' @return [`ggplot2::ggplot2`] object
#' @examples
#' data('scMethrix_data')
#' plot_coverage(scm = scMethrix_data)
#' @export
plot_coverage <- function(scm = NULL, type = c("histogram", "density"), pheno = NULL,
                          max_cov = 100, obs_lim = 1e+06, paletteID = "Dark Mint", show_legend = FALSE, verbose = TRUE,...) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  type <- .validateArg(type, plot_coverage)
  .validateType(pheno,c("string","null"))
  .validateType(max_cov,"integer")
  .validateType(obs_lim,"integer")
  .validateType(paletteID,"string")
  .validateType(show_legend,"boolean")
  
  Value <- Sample <- Pheno <- NULL
  
  colors_palette <- .getPalette(ncol(scm), paletteID = paletteID)
  
  #---- Function code ------------------------------------------------------
  if (matrixStats::product(dim(scm)) > obs_lim) {
    message("The dataset is bigger than the size limit. A random subset of the object will be used that contains ~",
            obs_lim, " observations.")
    n_rows <- trunc(obs_lim/ncol(scm))
    sel_rows <- sample(seq_len(nrow(scm)), size = n_rows,
                       replace = FALSE)
  } else {
    sel_rows <- seq_len(nrow(scm))
  }

  plot.data <- prepare_plot_data(scm = scm[sel_rows, ], assay = "counts", pheno = pheno, na.rm = T)
  #setnafill(plot.data,fill=0,cols="Value")
  
  plot.data <- plot.data[Value <= max_cov, ]

  # generate the plots
  if (is.null(pheno)) {
    if (type == "density") {
      p <- ggplot2::ggplot(plot.data, aes(Value, color = Sample)) +
        ggplot2::geom_density(alpha = 0.5, adjust = 3, lwd = 1.5, show.legend = show_legend, na.rm = T,
                              position = "identity") + ggplot2::theme_classic() +
        ggplot2::scale_color_manual(values = colors_palette) + 
        ggplot2::ylab("Density") + ggplot2::xlab("Coverage") +
        scale_x_sqrt(breaks=c(0,1,2,4,8))# + scale_y_sqrt(breaks=c(0,1,2,4,8))
      
    } else if (type == "histogram") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(Value, fill = Sample)) + 
        ggplot2::geom_histogram(alpha = 0.6, binwidth = 1, color = "black", show.legend = show_legend) + 
        ggplot2::theme_classic() +
        ggplot2::xlab("Coverage") + ggplot2::ylab("Density")
        ggplot2::scale_color_manual(values = colors_palette)
      # print(p)
    }
  } else {
    if (type == "density") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(Value, color = Pheno)) +
        ggplot2::geom_density(alpha = 0.6, adjust = 1.5, lwd = 1, show.legend = show_legend,
                              position = "identity") + ggplot2::theme_classic() + 
        ggplot2::xlab("Coverage") + ggplot2::ylab("Density") +
        ggplot2::scale_fill_manual(values = colors_palette)
      # print(p)
    } else if (type == "histogram") {
      p <- ggplot2::ggplot(plot.data, ggplot2::aes(Value, fill = Pheno)) +
        ggplot2::geom_histogram(alpha = 0.6, binwidth = 1, color = "black", show.legend = show_legend) + 
        ggplot2::theme_classic() + ggplot2::xlab("Coverage") + ggplot2::ylab("Density") +
        ggplot2::scale_fill_manual(values = colors_palette)
      # print(p)
    }
  }
  
  gc(verbose = FALSE)
  
  p <- p + ggplot2::theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12,
                                                                       colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),
            axis.title.y = element_blank(), legend.title = element_blank())
  
  return(p + scMethrix_theme(...))
}

#---- plot_sparsity ----------------------------------------------------------------------------------------------------
#' Sparsity of sample
#' @inheritParams plot_violin
#' @param type `string`; Choose between `Boxplot` or `Scatterplot`. Default = `Scatterplot`
#' @param by `string`; The variable for x-axis, `Sample` or `Chromosome.` If by `Sample`, and phenotype is not `NULL`, phenotype will be taken instead. Default = `Sample`
#' @param phenotype `string`; Col name of `colData(m)`. Will be used as a factor to color different groups
#' @param show_avg `boolean`; flag to show a dotted line representing the average value.
#' @return [`ggplot2::ggplot2`] object
#' @examples
#' data('scMethrix_data')
#' plot_sparsity(scm = scMethrix_data)
#' @export
plot_sparsity <- function(scm = NULL, assay = "score", type = c("Scatterplot", "Boxplot", "Jitterplot"), by = c("Sample","Chromosome"), phenotype = NULL, show_legend = FALSE, verbose = TRUE, show_avg = TRUE, paletteID = "Dark Mint",...) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  type <- .validateArg(type, plot_sparsity)
  by <- .validateArg(by, plot_sparsity)
  .validateAssay(scm,assay)
  .validateType(phenotype,c("string","null"))
  
  Sparsity <- variable <- NULL
  
  colors_palette <- .getPalette(ncol(scm), paletteID = paletteID)
  
  if (!is.null(phenotype) && type == "Scatterplot" && by == "Chromosome") {
    warning("Phenotype given for scatterplot when graphing by chromosome. Phenotype will be ignored.")
    phenotype = NULL
  }

  if (!is.null(phenotype) && by == "Sample") {
    if (!phenotype %in% colnames(colData(scm))) {
      stop("Please provide a valid phenotype annotation column. The column must exist within colData(scm).")
    } else {
      by = "Phenotype"
    }
  }
  
  Sample <- Phenotype <- Count <- Sites <- Chromosome <- . <- NULL
  
  #---- Function code ------------------------------------------------------
  chrs = rowRanges(scm)@seqnames
  end = cumsum(chrs@lengths)
  start = c(1, head(end, -1) + 1)
  chrs = data.frame(Chromosome = as.character(chrs@values), Start = start, End = end, Sites = end-start)
  
  stats <- getStats(scm,perChr=TRUE, perSample = F, phenotype = phenotype, stats="Count")
  stats <- merge(stats,chrs[,c("Chromosome","Sites")],by="Chromosome")
  #stats[,Sparsity := Count/Sites]
  
  if (!is.null(phenotype)) {
    pheno <- scm@colData[, phenotype,drop=FALSE]
    pheno$Sample = row.names(pheno)
    pheno <- as.data.table(pheno)
    #pheno[ , Members := .N, by = .(Group)]
    stats <- merge(stats,pheno,by="Sample")
    setnames(stats, phenotype, "Phenotype")
  } else {
    stats[,Phenotype := Sample]
  }
  
  stats <- stats[, lapply(.SD, as.character), by=.(Count,Sites)]
  x_lab <- "Sample"
  
  
  # Setup the graph data
  if (type == "Scatterplot") {
    
    sparsity <- stats[, .(Sparsity = sum(Count)/sum(Sites)*100), by = by]
    gg_data <- ggplot2::ggplot(sparsity, aes_string(x = by, y = "Sparsity", color = by))
    gg_elem <- ggplot2::geom_point(show.legend = show_legend)
    
  } else {
    
    if (!is.null(phenotype)) {
      sparsity <- stats[, .(Sparsity = sum(Count)/sum(Sites)*100), by = list(Phenotype,Chromosome)]
    } else {
      sparsity <- stats[, .(Sparsity = sum(Count)/sum(Sites)*100), by = list(Sample,Chromosome)]
    }
    
    if (type == "Boxplot") {
      gg_elem <- ggplot2::geom_boxplot(show.legend = show_legend)
      gg_color <- by
      gg_data <- ggplot2::ggplot(sparsity, aes_string(x=by, y="Sparsity", fill = gg_color))
    } else {
      dodge_dist = 2#length(ggplot_build(gg_data)$plot$scales$scales[[3]]$get_labels())
      gg_elem <- ggplot2::geom_point(show.legend = show_legend, 
                                     position = position_dodge(width = 1/dodge_dist))
      if (by == "Phenotype" || by == "Sample") gg_color <- "Chromosome"
      else if (!is.null(phenotype)) gg_color <- "Phenotype"
      else gg_color <- "Sample"
      gg_data <- ggplot2::ggplot(sparsity, aes_string(x=by, y="Sparsity", color = gg_color))
    }
  }

  p <- gg_data + gg_elem + 
    ggplot2::scale_color_manual(values = .getPalette(ncol(scm))) + #,name = legend_lab) + 
    ylab("Sparsity (%)")# + xlab(x_lab) 
  
  if (show_avg) {
    avg.spars = mean(sparsity$Sparsity)
    p <- p + geom_hline(yintercept=avg.spars, linetype="dashed",
                       color = "black", size=1)+
            scale_y_continuous(
              sec.axis = dup_axis(
                breaks = avg.spars,
                labels = parse(text="bar(x)"),
                name = NULL
              )
            )
  }
  
  return (p + scMethrix_theme(...))
}

# plotChrStats <- function(scm, assay = "score", stat = c("Mean", "Median", "Count", "Proportion"), perSample = FALSE, 
#                          phenotype = NULL, ignoreChrs = NULL, ignoreSamples = NULL, verbose = TRUE, show_legend = FALSE,...) {
#   
#   getStat <- if(stat == "Mean" || stat == "Median") c(stat,"SD") else stat
#   plotData <- getStats(scm, assay = assay, stats = getStat, perSample = perSample, perChr = TRUE, 
#                        ignoreChrs = ignoreChrs, ignoreSamples = ignoreSamples)
#   
#   .plotStats(plotData = plotData, by = "Chromosome", collapsed = !perSample)
#   
# }
# 
# plotSampleStats <- function(scm, assay = "score", stat = c("Mean", "Median", "Count", "Proportion"), perChr = FALSE, 
#                             phenotype = NULL, ignoreChrs = NULL, ignoreSamples = NULL, verbose = TRUE, show_legend = FALSE,...) {
#   
#   getStat <- if(stat == "Mean" || stat == "Median") c(stat,"SD") else stat
#   plotData <- getStats(scm, assay = assay, stats = getStat, perSample = TRUE, perChr = perChr, 
#                        ignoreChrs = ignoreChrs, ignoreSamples = ignoreSamples)
#   
#   .plotStats(plotData = plotData, by = "Sample", collapsed = !perChr)
# }


#---- plotStats -------------------------------------------------------------------------------------------------------
#' Plot descriptive statistics results from [getStats()]
#' @details Plot descriptive statistics
#' @inheritParams plot_violin
#' @param stat `string`; Can plot `Mean`, `Median`, `Count`, or `Proportion`. `Count` is number of non-NA assay values, whereas `Proportion` is the proportion of non-NA assay values. Default = `Mean`.
#' @param type `string`; Choose between `Boxplot` or `Scatterplot`. Only applies when `collapse = FALSE`. Default = `Scatterplot`.
#' @param by `string`; Can plot via `Sample` or `Chromosome`. Default = `Sample`.
#' @param collapse `boolean`; Collapse by sample or chromosome. Will collapse by the opposite of `by` value. Default = `TRUE`.
#' @param phenotype `string`; Group samples by `phenotype`. The `phenotype` must have a corresponding column in `colData()`. Default = `NULL`.
#' @param ignoreChrs `list(string)`; Chromosomes to ignore. Default = `NULL`.
#' @param ignoreSamples `list(string)`; Samples to ignore. Default = `NULL`.
#' @param nCol `integer`; Number of columns. Passed to [ggplot2::facet_wrap()].
#' @param nRow `integer`; Number of rows. Passed to [ggplot2::facet_wrap()].
#' @return [`ggplot2::ggplot2`] object
#' @seealso [getStats()]
#' @examples
#' data('scMethrix_data')
#' plotStats(scMethrix_data)
#' @export
#'
plotStats <- function(scm, assay = "score", stat = c("Mean", "Median", "Count", "Proportion"), 
                      type = c("Scatterplot", "Boxplot"), by = c("Sample", "Chromosome"), collapse = TRUE, 
                      phenotype = NULL, ignoreChrs = NULL, ignoreSamples = NULL, nCol = NULL, nRow = NULL, 
                      verbose = TRUE, show_legend = FALSE, paletteID = "Dark Mint",...) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  stat <- .validateArg(stat,plotStats)
  type <- .validateArg(type,plotStats)
  by <- .validateArg(by,plotStats)
  .validateType(phenotype,c("string","null"))
  .validateType(ignoreChrs,c("string","null"))
  .validateType(ignoreSamples,c("string","null"))
  .validateType(nCol,c("integer","null"))
  .validateType(nRow,c("integer","null"))
  .validateType(verbose,"boolean")
  
  if (!is.null(phenotype) && !phenotype %in% colnames(colData(scm))) 
    stop("Error in plotting. No column named '",phenotype,"' is present in colData(). Must be one of: ",
         .pasteList(colnames(colData(scm))))

  if (any(!ignoreSamples %in% sampleNames(scm)))
    warning("Ignored samples are not present in the data: ",
            .pasteList(ignoreSamples[!ignoreSamples %in% sampleNames(scm)]))
  
  if (any(!ignoreChrs %in% levels(seqnames(scm))))
    warning("Ignored chromosomes are not present in the data: ",
            .pasteList(ignoreChrs[!ignoreChrs %in% levels(seqnames(scm))]))
  
  Chromosome <- . <- Sample <- Proportion <- SD <- SDlow <- SDhigh <- Value <- NULL

  #- Function code --------------------------------------------------------------------------
  
  yTitle <- tools::toTitleCase(stat)
  
  colorPalette <- .getPalette(ncol(scm), paletteID = paletteID)
  
  getStat <- if (stat == "Mean" || stat == "Median") c(stat,"SD") else "Count"
  perChr <- if (by == "Chromosome") TRUE else !collapse
  perSample <- if (by == "Sample") TRUE else !collapse
  group <- if (by == "Chromosome") "Sample" else "Chromosome" 
  
  plotData <- getStats(scm, assay = assay, stats = getStat, perSample = perSample, perChr = perChr, 
                       ignoreChrs = ignoreChrs, ignoreSamples = ignoreSamples)
  setDT(plotData)
 
  if (stat == "Proportion") {
    setnames(plotData, "Count", "Proportion")
    plotData[, `:=`(Proportion, Proportion/nrow(scm))]
  }

  # plot_dat$Chr <- gsub("^.{0,3}", "", plot_dat$Chr)
  # plot_dat$Chr <- factor(plot_dat$Chr, levels = unique(plot_dat$Chr))
   
  setnames(plotData, stat, "Value")

  figure <- ggplot(data = plotData, aes_string(x = by, y = "Value", label = group)) +
    ggplot2::ylab(stat) + 
    ggplot2::xlab(by) +
    ggplot2::theme_minimal(base_size = 12) + 
    ggplot2::theme(axis.text.x = element_text(hjust = 1, size = 10, colour = "black"),
                   axis.text.y = element_text(size = 10, colour = "black")) +
   scMethrix_theme()

  if (collapse) {
    figure <- figure + geom_point()
    if (stat == "Mean") {
      plotData[, `:=`(SD, as.numeric(as.character(SD)))]
      plotData[, `:=`(SDlow,  max(Value - SD,0))]
      plotData[, `:=`(SDhigh, Value + SD)]
      figure <- figure + ggplot2::geom_errorbar(aes(ymin = SDlow, ymax = SDhigh), col = "gray25")
    }
  } else {
    nPoints <- length(unique(plotData[[group]])) #TODO: this should check the variable with fewest groups
    if (nPoints >= 5 && type == "Boxplot") {
      figure <- figure + ggplot2::geom_boxplot()
    } else {
      figure <- figure + ggplot2::geom_jitter(width = .calcJitter(nPoints))
    }
  }
    
  #plotly::ggplotly(p = figure)
  
  return(figure  + scMethrix_theme(...))
}

#---- plot_dim_red -----------------------------------------------------------------------------------------------------
#' Plot dimensionality reduction
#' @inheritParams generic_scMethrix_function
#' @inheritParams plot_violin
#' @param dim_red `string`; name of a dimensionality reduction inside an [`scMethrix`] object. Should be a matrix of two columns representing
#' the X and Y coordinates of the dim. red., with each row being a separate sample
#' @param axis_labels `list(string)`; A named list of `X` and `Y` strings for labels, or `NULL` if no labels are desired
#' @param color_anno `string`; Column name of `colData(scm)`. Default = `NULL.` Will be used as a factor to color different groups.
#' @param shape_anno string; Column name of `colData(scm)`. Default = `NULL.` Will be used as a factor to shape different groups.
#' @param show_dp_labels boolean; Flag to show the labels for dots. Default = `FALSE`
#' @return [`ggplot2::ggplot2`] object
#' @importFrom graphics par mtext lines axis legend title
#' @export
# plot_dim_red <- function(scm, dim_red, palette = "Paired", color_anno = NULL, shape_anno = NULL, legend_anno = NULL, axis_labels = NULL, show_dp_labels = FALSE, verbose = TRUE) {
#   
#   #---- Input validation ---------------------------------------------------
#   X <- Y <- Color <- Shape <- color <- shape <- shapes  <- colors <- Sample <- row_names <- NULL
#   
#   .validateExp(scm)
#   .validateType(dim_red,"string")
#   if (!(dim_red %in% reducedDimNames(scm))) stop("Invalid dim_red specified. '",dim_red,"' does not exist in the experiment.")
#   .validateType(color_anno,c("string","null"))
#   .validateType(shape_anno,c("string","null"))
#   .validateType(axis_labels,c("string","null"))
#   .validateType(show_dp_labels,"boolean")
#   
#   dim_red <- reducedDim(scm,type=dim_red)
#   
#   if (ncol(dim_red) != 2) {
#     warning("More than two columns in the dimentionality reduction. Only the first two will be used")
#     dim_red <- dim_red[,1:2]
#   }
# 
#   dim_red = as.data.frame(dim_red)
#   colnames(dim_red) <- c("X", "Y")
#   dim_red <- merge(dim_red,colData(scm),by="row.names")
#   names(dim_red)[names(dim_red) == 'Row.names'] <- 'Sample'
#   dim_red$Sample <- as.character(dim_red$Sample)
#   
#   dim_red <- as.data.frame(dim_red)
#   #dim_red[, shape_anno] <- as.factor(dim_red[, shape_anno])
#   
#   #---- Function code ------------------------------------------------------
#   
#   # if (!is.null(color_anno)) {
#   #   if (color_anno  %in% colnames(colData(scm))) {
#   #     dim_red$Color <- as.factor(unlist(as.data.table(colData(scm))[,color_anno, with=FALSE]))
#   #     colors <- scale_color_manual(values= .getPalette(length(unique(dim_red$Color)),palette = palette))
#   #   } else {
#   #     stop(paste0(color_anno, " not found in provided scMethrix object"))
#   #   }
#   # }
#   # 
#   # if (!is.null(shape_anno)) {
#   #   if (shape_anno %in% colnames(colData(scm))) {
#   #     dim_red$Shape <- as.factor(unlist(as.data.table(colData(scm))[,shape_anno, with=FALSE])) 
#   #     shapes <- scale_shape_manual(values = .getShapes(length(unique(dim_red$Shape))))
#   #   } else {
#   #     stop(paste0(shape_anno, " not found in provided scMethrix object"))
#   #   }
#   # }  
#   
#   if (is.null(axis_labels)) {
#     axis_labels = list(X="",Y="")
#   }
#   
#   
#   Cell_order <- dim_red$Order
#   names(Cell_order) <- dim_red$Cell
#   Cell_order <- Cell_order[!duplicated(names(Cell_order))]
#   Cell_order <- names(sort(Cell_order))
#   dim_red$Cell <- factor(dim_red$Cell, levels = Cell_order)
#   
#   Cell_color <- dim_red$Color
#   names(Cell_color) <- dim_red$Cell
#   Cell_color <- Cell_color[!duplicated(names(Cell_color))]
#   Cell_color <- Cell_color[order(factor(names(Cell_color), levels=Cell_order))]
#   
#   Cell_shape <- dim_red$Shape
#   names(Cell_shape) <- dim_red$Cell
#   Cell_shape <- Cell_shape[!duplicated(names(Cell_shape))]
#   Cell_shape <- Cell_shape[order(factor(names(Cell_shape), levels=Cell_order))]
# 
#   browser()
#   
#   if (all(c("Color", "Shape") %in% colnames(dim_red))) {
#   #  dimred_gg <- 
#      dimred_gg <-  ggplot2::ggplot(data = dim_red, aes(x = X, y = Y)) + 
#       geom_point(aes(fill = Cell, shape = Cell),size = 4, stroke = .2, alpha = .5) +
#       scale_shape_manual(values = Cell_shape) + 
#       scale_fill_manual(values = Cell_color) 
#       
#      # geom_mark_hull(aes(fill=Cell),concavity = 1)
#       
#       
#     #   
#     #   scale_shape_manual(values = Cell_shape, labels = names(Cell_shape)) +
#     # scale_fill_manual(values = Cell_color, labels = names(Cell_color))
#     #   
#     #   
#     #   
#     #   scale_discrete_manual_ext(c("group", "color","shape"), values = list(color = c("yellow", "darkred"),linetype = c("solid", "dashed")), name = "Legend name")
#     #   
#     #   
#     #   
#     # scale_discrete_manual(aes(group = Cell, color = Color, shape = Shape))
#     # 
#     # 
#     # 
#     #    +
#     #     scale_group_discrete("Model 1") +
#     # 
#     #   
#     #   
#     #   
#     #   
#     # labs(color  = "Cell type", group = "Cell type", shape = "Cell type")
#       
#       
#   } else if ("Color" %in% colnames(dim_red)) {
#     dimred_gg <- ggplot2::ggplot(data = dim_red, aes(x = X, y = Y, color = Color,
#                                                      label = Sample))  + labs(color = color_anno)
#   } else if ("Shape" %in% colnames(dim_red)) {
#     dimred_gg <- ggplot2::ggplot(data = dim_red, aes(x = X, y = Y, shape = Shape,
#                                                      label = Sample)) + labs(shape = shape_anno)
#   } else {
#     dimred_gg <- ggplot2::ggplot(data = as.data.frame(dim_red), aes(x = X, y = Y,
#                                                                     label = Sample))
#   }
# 
#   dimred_gg <- dimred_gg  + ggplot2::theme_classic(base_size = 12) + 
#     ggplot2::xlab(axis_labels$X) + ggplot2::ylab(axis_labels$Y) + 
#     ggplot2::theme(axis.text.x = element_text(colour = "black", size = 12),
#                    axis.text.y = element_text(colour = "black", size = 12)) 
#   #+ geom_point(size = 4, stroke = .2, alpha = .5)
#   
#   if (show_dp_labels) dimred_gg <- dimred_gg + ggplot2::geom_label(size = 4) 
#   
#   return(dimred_gg)
#   
# }
plot_dim_red <- function(scm, dim_red, paletteID = "Dark Mint", color_anno = NULL, shape_anno = NULL, axis_labels = NULL, show_dp_labels = FALSE, verbose = TRUE,...) {

  #---- Input validation ---------------------------------------------------
  X <- Y <- Color <- Shape <- color <- shape <- shapes  <- colors <- Sample <- row_names <- NULL

  .validateExp(scm)
  .validateType(dim_red,"string")
  if (!(dim_red %in% reducedDimNames(scm))) stop("Invalid dim_red specified. '",dim_red,"' does not exist in the experiment.")
  .validateType(color_anno,c("string","null"))
  .validateType(shape_anno,c("string","null"))
  .validateType(axis_labels,c("string","null"))
  .validateType(show_dp_labels,"boolean")

  dim_red <- reducedDim(scm,type=dim_red)

  if (ncol(dim_red) != 2) {
    warning("More than two columns in the dimentionality reduction. Only the first two will be used")
    dim_red <- dim_red[,1:2]
  }

  dim_red = as.data.frame(dim_red)
  colnames(dim_red) <- c("X", "Y")
  dim_red$Sample = rownames(dim_red)

  #---- Function code ------------------------------------------------------

  if (!is.null(color_anno)) {
    if (color_anno  %in% colnames(colData(scm))) {
      dim_red$Color <- as.factor(unlist(as.data.table(colData(scm))[,color_anno, with=FALSE]))
      colors <- scale_color_manual(values= .getPalette(length(unique(dim_red$Color)),paletteID = paletteID))
    } else {
      stop(paste0(color_anno, " not found in provided scMethrix object"))
    }
  }

  if (!is.null(shape_anno)) {
    if (shape_anno %in% colnames(colData(scm))) {
      dim_red$Shape <- as.factor(unlist(as.data.table(colData(scm))[,shape_anno, with=FALSE]))
      shapes <- scale_shape_manual(values = .getShapes(length(unique(dim_red$Shape))))
    } else {
      stop(paste0(shape_anno, " not found in provided scMethrix object"))
    }
  }

  if (is.null(axis_labels)) {
    axis_labels = list(X="",Y="")
  }

  if (all(c("Color", "Shape") %in% colnames(dim_red))) {
    dimred_gg <- ggplot2::ggplot(data = dim_red, aes(x = X, y = Y, color = Color,
                                            shape = Shape, label = Sample)) +
      labs(color = color_anno, shape = shape_anno)
  } else if ("Color" %in% colnames(dim_red)) {
    dimred_gg <- ggplot2::ggplot(data = dim_red, aes(x = X, y = Y, color = Color,
                                            label = Sample))  + labs(color = color_anno)
  } else if ("Shape" %in% colnames(dim_red)) {
    dimred_gg <- ggplot2::ggplot(data = dim_red, aes(x = X, y = Y, shape = Shape,
                                            label = Sample)) + labs(shape = shape_anno)
  } else {
    dimred_gg <- ggplot2::ggplot(data = as.data.frame(dim_red), aes(x = X, y = Y,
                                                           label = Sample))
  }

  dimred_gg <- dimred_gg  + ggplot2::theme_classic(base_size = 12) +
    ggplot2::xlab(axis_labels$X) + ggplot2::ylab(axis_labels$Y) +
    ggplot2::theme(axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 12)) + geom_point(size = 4, stroke = .2, alpha = .5)

  if (show_dp_labels) dimred_gg <- dimred_gg + ggplot2::geom_label(size = 4)

  dimred_gg <- dimred_gg + colors + shapes + scMethrix_theme(...)

  return(dimred_gg)

}

#------------------------------------------------------------------------------------------------------------
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

#---- benchmark_imputation ---------------------------------------------------------------------------------------------
#' Evaluates imputations methods by NRMSE or AUC
#' @details Does stuff
#' @param sparse_prop `numeric`; A sparsity proportion between 0 and 1. E.g. 0.1 replaces 10% of the matrix with NA
#' @param imp_methods `closure`; The imputation methods to compare.
#' @param iterations `integer`; Number of iterations to test
#' @param type `character`; descriptive statistic. Can be either `AUC` or `RMSE.` Default `RMSE`
#' @inheritParams generic_scMethrix_function
#' @return [`ggplot2::ggplot2`]; The graph showing the NRMSE for each imputation method at each sparsity
#' @examples
#' data('scMethrix_data')
#' \dontrun{
#' scMethrix_data <- impute_regions(scMethrix_data, new_assay="impute",type="RF")
#' benchmark_imputation(scMethrix_data, assay="impute", sparse_prop = c(0.1,0.5,0.85))
#' }
#' @export
benchmark_imputation <- function(scm = NULL, assay = "score", sparse_prop = seq(0.1, 0.9, 0.1), iterations = 3,
                                 imp_methods = c(iPCA = function(...) impute_regions(type="iPCA",...), 
                                                 RF = function(...) impute_regions(type="RF",...), 
                                                 kNN = function(...) impute_regions(type="kNN",...)),
                                 type = "RMSE") {

  #---- Input validation ---------------------------------------------------
  .validatePackageInstall("Metrics")
  
  . <- results <- Sparsity <- NRMSE <- Imputation <- AUC <- NULL
  
  #---- Function code ------------------------------------------------------
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
    ggplot2::scale_x_continuous(breaks=seq(0,1,.05))  
    #scMethrix_theme() 
}

#---- scMethrix_theme --------------------------------------------------------------------------------------------------
#' Theme for ggplot
#' @param base_size `integer`; Size of text
#' @param base_family `string`; Family of text
#' @param ... Additional arguments for [ggplot2::theme()]
#' @return [`ggplot2::ggplot2`] element; data for the ggplot theme
#' @export
scMethrix_theme <- function(base_size = 12, base_family = "",...) {

  update_geom_defaults("line", list(size = 1.2))
  update_geom_defaults("point", list(size = 3))

  theme_bw(base_size = base_size, base_family = base_family) %+replace%
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
      axis.text.x = element_text(angle = 90),
      axis.line =         element_line(colour = '#969696', size = 1, #TODO: Prefer only bottom line
                                       linetype = 1, lineend = "butt"),
      panel.grid.major =  element_line(colour = '#DADADA', size = 0.75,
                                       linetype = 1, lineend = "butt"),
      panel.grid.minor =  element_blank(),
      plot.background =   element_blank(),
      panel.background =  element_rect(),
      legend.key =        element_rect(colour = '#DADADA'),...,
      complete = TRUE
    )
}


#---- .calcJitter ------------------------------------------------------------------------------------------------------
#' Calculates jitter based on number of data points
#' @description Jittering of few points can look bad, depending on the RNG, as some groups will have points on terminal edges of the space, while others will be nearly in line. This limits the spread of the jitter to improve visualization.
#' @param n integer; Number of data points
#' @param max numeric; Maximum jitter
#'
#' @return a number between 0 and 1
#' @nord
#' @examples
#' .calcJitter(5)
.calcJitter <- function(n, max = 0.8) {
  max <- max(min(1,max),0)
  return(max(min(max,n*0.05),0))
}