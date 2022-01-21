test_that("prepare_plot_data", {
  expect_error(get_region_summary("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(prepare_plot_data(scm, n_cpgs = "not an int"),msg.validateType)
    
    d <- prepare_plot_data(scm,na.rm=F)
    expect_equal(dim(d),c(nrow(scm)*ncol(scm),3))
    expect_equal(colnames(d),c("Sample","Value","Pheno"))
    
    n_NAs <- nrow(d[is.na(Value)]) 
    d <- prepare_plot_data(scm,)
    expect_equal(dim(d),c(nrow(scm)*ncol(scm)-n_NAs,3))
    expect_equal(colnames(d),c("Sample","Value","Pheno"))
    
    n_cpgs = 50
    d <- prepare_plot_data(scm,n_cpgs = n_cpgs,na.rm=F)
    expect_equal(dim(d),c(n_cpgs*ncol(scm),3))
    expect_equal(colnames(d),c("Sample","Value","Pheno"))
    
    n_NAs <- nrow(d[is.na(Value)]) 
    d <- prepare_plot_data(scm,n_cpgs = n_cpgs)
    expect_equal(dim(d),c(n_cpgs*ncol(scm)-n_NAs,3))
    expect_equal(colnames(d),c("Sample","Value","Pheno"))
    
  }))
})

invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  
  test_that("plot_violin", {graph_test_helper(scm,plot_violin)})

  test_that("plot_density", {graph_test_helper(scm,plot_density)})

  test_that("plot_coverage", {
    graph_test_helper(scm, plot_coverage, type="histogram")
    graph_test_helper(scm, plot_coverage, type="density")
    graph_test_helper(scm, plot_coverage, type="histogram", pheno="Group")
    graph_test_helper(scm, plot_coverage, type="density",   pheno="Group")
  })
  
  test_that("plot_sparsity", {

      samp <- sampleNames(scm)
      chr <- levels(rowRanges(scm)@seqnames)
      pheno <- unique(colData(scm)$Group)
    
      plot = plot_sparsity(scm,type = "Scatterplot", by = "Sample")
      graph_test_helper2(plot, expected_x = samp)
      plot = plot_sparsity(scm,type = "Scatterplot", phenotype = "Group")
      graph_test_helper2(plot, expected_x = pheno)
      plot = plot_sparsity(scm,type = "Scatterplot", by = "Chromosome")
      graph_test_helper2(plot, expected_x = chr)
      
      plot = plot_sparsity(scm,type = "Boxplot", by = "Sample")
      graph_test_helper2(plot, expected_x = samp)
      plot = plot_sparsity(scm,type = "Boxplot", phenotype = "Group")
      graph_test_helper2(plot, expected_x = pheno)
      plot = plot_sparsity(scm,type = "Boxplot", by = "Chromosome")
      graph_test_helper2(plot, expected_x = chr)
      
      plot = plot_sparsity(scm,type = "Jitterplot", by = "Sample")
      graph_test_helper2(plot, expected_x = samp, expected_group = chr)
      plot = plot_sparsity(scm,type = "Jitterplot", phenotype = "Group")
      graph_test_helper2(plot, expected_x = pheno, expected_group = chr)
      
      plot = plot_sparsity(scm,type = "Jitterplot", by = "Chromosome")
      graph_test_helper2(plot, expected_x = chr, expected_group = samp)
      plot = plot_sparsity(scm,type = "Jitterplot", by = "Chromosome", phenotype = "Group",show_legend = T)
      graph_test_helper2(plot, expected_x = chr, expected_group = pheno)
      
  })
  
  test_that("plot_stats", {
    graph_test_helper(scm, plot_stats, per_chr = F, indiv_chr = F)
    graph_test_helper(scm, plot_stats, per_chr = T, indiv_chr = T)
  })
  
  test_that("plot_dim_red", {
    invisible(lapply(list("PCA","tSNE","UMAP"), function(type) {
    
      scm.dimred <- suppressWarnings(impute_regions(scm))
      scm.dimred <- dim_red_scMethrix(scm.dimred, assay="impute", type = type, verbose = FALSE) 
      
      dim_red_graph_test_helper(scm.dimred, plot_dim_red, dim_red=type)
      dim_red_graph_test_helper(scm.dimred, plot_dim_red, dim_red=type, color_anno="Group", shape_anno="Group")
    }))
  })
}))
