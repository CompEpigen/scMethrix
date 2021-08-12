test_that("reduce_cpgs", {
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_error(reduce_cpgs(scm="not scMethrix"),"A valid scMethrix object needs to be supplied")
    expect_error(reduce_cpgs(scm,assay="not an assay"),"Assay does not exist in the object")
    expect_error(reduce_cpgs(scm,var="not a var"),"'arg' should be one of")
    expect_error(reduce_cpgs(scm,top_var = 0,var="top"),"Zero loci available post NA removal")
    
    expect_equal(reduce_cpgs(scm,top_var = NULL),score(scm))
    
    cpgs = 10
    expect_equal(dim(reduce_cpgs(scm,top_var = cpgs,var="rand")),c(cpgs,ncol(scm)))
    expect_equal(dim(reduce_cpgs(scm,top_var = cpgs,var="top")),c(cpgs,ncol(scm)))

  }))
})

test_that("dim_red_scMethrix", {

  invisible(lapply(list(scm.mem,scm.h5), function(scm) {

    #PCA
    n_pc = 2
    pca <- dim_red_scMethrix(scm, n_pc = n_pc,type="PCA")
    expect_equal(reducedDimNames(pca),"PCA")
    expect_equal(dim(reducedDim(pca)),c(ncol(scm),n_pc))
    expect_equal(length(pca@metadata$PCA_vars),n_pc)
    
    #UMAP
    umap <- dim_red_scMethrix(scm,type="UMAP")
    expect_equal(reducedDimNames(umap),"UMAP")
    expect_equal(dim(reducedDim(umap)),c(ncol(scm),2))
    
    #tSNE
    tsne <- dim_red_scMethrix(scm,type="tSNE")
    expect_equal(reducedDimNames(tsne),"tSNE")
    expect_equal(dim(reducedDim(tsne)),c(ncol(scm),2))
  }))
})