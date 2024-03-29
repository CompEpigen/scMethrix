test_that("prepare_plot_data", {
  expect_error(get_region_summary("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(prepare_plot_data(scm, n_cpgs = "not an int"),msg.validateType)
    expect_error(prepare_plot_data(scm, regions = "not a range"),msg.validateType)
    
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
    graph_test_helper(scm, plot_sparsity, indiv_samples = F, type="box")
    graph_test_helper(scm, plot_sparsity, indiv_samples = F, type="scatter")
    graph_test_helper(scm, plot_sparsity, indiv_samples = F, type="box",     pheno="Group")
    graph_test_helper(scm, plot_sparsity, indiv_samples = F, type="scatter", pheno="Group", )
  })
  
  test_that("plot_stats", {
    graph_test_helper(scm, plot_stats, per_chr = F, indiv_chr = F)
    graph_test_helper(scm, plot_stats, per_chr = T, indiv_chr = T)
  })
  
  test_that("plot_dim_red", {
    invisible(lapply(list("PCA","tSNE","UMAP"), function(type) {
    
      scm.dimred <- dim_red_scMethrix(scm,type = type, verbose = FALSE) 
      
      dim_red_graph_test_helper(scm.dimred, plot_dim_red, dim_red=type)
      dim_red_graph_test_helper(scm.dimred, plot_dim_red, dim_red=type, color_anno="Group", shape_anno="Group")
    }))
  })
}))
