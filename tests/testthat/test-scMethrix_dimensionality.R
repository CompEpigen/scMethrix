test_that("reduce_cpgs", {
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_error(reduce_cpgs(scm="not scMethrix"))
    expect_error(reduce_cpgs(scm,assay="not an assay"))
    expect_error(reduce_cpgs(scm,var="not a var"))
    expect_error(reduce_cpgs(scm,top_var = 0,var="top"))
    
    expect_equal(reduce_cpgs(scm,top_var = NULL),score(scm))
    
    cpgs = 10
    expect_equal(dim(reduce_cpgs(scm,top_var = cpgs,var="rand")),c(cpgs,ncol(scm)))
    expect_equal(dim(reduce_cpgs(scm,top_var = cpgs,var="top")),c(cpgs,ncol(scm)))

  }))
})

test_that("pca_scMethrix", {
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {

    n_pc = 2
    pca <- pca_scMethrix(scm, n_pc = n_pc)
    expect_equal(reducedDimNames(pca),"PCA")
    expect_equal(dim(reducedDim(pca)),c(ncol(scm),n_pc))
    expect_equal(length(pca@metadata$PCA_vars),n_pc)
  }))
})

test_that("umap_scMethrix", {
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    umap <- umap_scMethrix(scm)
    expect_equal(reducedDimNames(umap),"UMAP")
    expect_equal(dim(reducedDim(umap)),c(ncol(scm),2))
    
  }))
})

test_that("tsne_scMethrix", {
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    tsne <- tsne_scMethrix(scm)
    expect_equal(reducedDimNames(tsne),"tSNE")
    expect_equal(dim(reducedDim(tsne)),c(ncol(scm),2))
    
  }))
})

