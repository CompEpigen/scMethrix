test_that("prepare_plot_data", {
  expect_error(prepare_plot_data("not scMethrix"))
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(prepare_plot_data(scm, n_cpgs = "not an int"))
    expect_error(prepare_plot_data(scm, ranges = "not a range"))
    
    d <- prepare_plot_data(scm)
    expect_equivalent(dim(d),c(nrow(scm)*ncol(scm),2))
    expect_equivalent(colnames(d),c("variable","Meth"))
    
    n_cpgs = 50
    d <- prepare_plot_data(scm,n_cpgs = n_cpgs)
    expect_equivalent(dim(d),c(n_cpgs*ncol(scm),2))
    
  }))
})

test_that("get_palette", {
  expect_error(get_palette(0,"RdYlGn"))
  expect_error(get_palette(1,"not a palette"))
  
  colors = 5
  expect_length(get_palette(colors,"RdYlGn"),colors)
})

test_that("plot_violin", {
  expect_error(plot_violin("not scMethrix"))
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot = plot_violin(scm)
    expect_true("ggplot" %in% class(plot))
    expect_equivalent(levels(plot$data$variable),colnames(scm))
  }))
})

test_that("plot_density", {
  expect_error(plot_density("not scMethrix"))
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot = plot_density(scm)
    expect_true("ggplot" %in% class(plot))
    expect_equivalent(levels(plot$data$variable),colnames(scm))
  }))
})

test_that("plot_coverage", {
  expect_error(plot_coverage("not scMethrix"))
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot = plot_coverage(scm)
    expect_true("ggplot" %in% class(plot))
    expect_equivalent(levels(plot$data$variable),colnames(scm))
  }))
})

test_that("plot_stats", {
  expect_error(plot_stats("not scMethrix"))
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot = plot_stats(get_stats(scm))
    expect_true("ggplot" %in% class(plot))
    expect_equivalent(sort(unique(plot$data$Sample_Name)),sort(colnames(scm)))
    expect_equivalent(sort(unique(plot$data$Chromosome)),sort(levels(seqnames(scm))))
  }))
})

test_that("plot_pca", {
  expect_error(plot_pca("not scMethrix"))
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot <- plot_pca(pca_scMethrix(scm))
    expect_true("ggplot" %in% class(plot))
    expect_equivalent(plot$data$row_names,colnames(scm))
  }))
})

test_that("plot_tsne", {
  expect_error(plot_tsne("not scMethrix"))
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot <- plot_tsne(tsne_scMethrix(scm))
    expect_true("ggplot" %in% class(plot))
    expect_equivalent(plot$data$row_names,colnames(scm))
  }))
})

test_that("plot_umap", {
  expect_error(plot_umap("not scMethrix"))
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot <- plot_umap(umap_scMethrix(scm))
    expect_true("ggplot" %in% class(plot))
    expect_equivalent(plot$data$row_names,colnames(scm))
  }))
})
