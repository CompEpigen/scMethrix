test_that("prepare_plot_data", {
  expect_error(prepare_plot_data("not scMethrix") ,"A valid scMethrix object needs to be supplied")
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(prepare_plot_data(scm, n_cpgs = "not an int"),"n_cpgs must be numeric.")
    expect_error(prepare_plot_data(scm, ranges = "not a range"),"Invalid input class for regions. Must be a GRanges or data.frame-like")
    
    d <- prepare_plot_data(scm)
    expect_equal(dim(d),c(nrow(scm)*ncol(scm),2))
    expect_equal(colnames(d),c("variable","Meth"))
    
    n_cpgs = 50
    d <- prepare_plot_data(scm,n_cpgs = n_cpgs)
    expect_equal(dim(d),c(n_cpgs*ncol(scm),2))
    
  }))
})

test_that("get_palette", {
  expect_error(get_palette(0,"RdYlGn"),"Zero colors present in the palette")
  expect_error(get_palette(1,"not a palette"),"Please provide a valid RColorBrewer palettte.")
  
  colors = 5
  expect_length(get_palette(colors,"RdYlGn"),colors)
})

test_that("plot_violin", {
  graph_test_helper(plot_violin)
})

test_that("plot_density", {
  graph_test_helper(plot_density)
})

test_that("plot_coverage", {
  graph_test_helper(plot_coverage)
})

test_that("plot_sparsity", {
  graph_test_helper(plot_sparsity,indiv_samples = FALSE)
})

test_that("plot_stats", {
  expect_error(plot_pca("not scMethrix"),"A valid scMethrix object needs to be supplied")
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot = plot_stats(get_stats(scm))
    expect_true("ggplot" %in% class(plot))
    expect_equal(sort(unique(plot$data$Sample_Name)),sort(colnames(scm)))
    expect_equal(sort(unique(plot$data$Chromosome)),sort(levels(seqnames(scm))))
  }))
})

test_that("plot_pca", {
  expect_error(plot_pca("not scMethrix"),"A valid scMethrix object needs to be supplied")
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot <- plot_pca(pca_scMethrix(scm))
    expect_true("ggplot" %in% class(plot))
    expect_equal(plot$data$row_names,colnames(scm))
  }))
})

test_that("plot_tsne", {
  expect_error(plot_tsne("not scMethrix"),"A valid scMethrix object needs to be supplied")
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot <- plot_tsne(tsne_scMethrix(scm))
    expect_true("ggplot" %in% class(plot))
    expect_equal(plot$data$row_names,colnames(scm))
  }))
})

test_that("plot_umap", {
  expect_error(plot_umap("not scMethrix"),"A valid scMethrix object needs to be supplied")
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot <- plot_umap(umap_scMethrix(scm))
    expect_true("ggplot" %in% class(plot))
    expect_equal(plot$data$row_names,colnames(scm))
  }))
})