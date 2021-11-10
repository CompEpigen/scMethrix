test_that("reduce_scMethrix", {
  
  expect_error(reduce_scMethrix("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_error(reduce_scMethrix(scm="not scMethrix"),msg.validateExp)
    expect_error(reduce_scMethrix(scm,assay="not an assay"),msg.validateAssay)
    expect_error(reduce_scMethrix(scm,var="not a var"),msg.validateArg)
    expect_error(reduce_scMethrix(scm,n_cpg = "not an int"),msg.validateType)
    expect_error(reduce_scMethrix(scm,n_cpg = 0,var="top"),"Zero loci available post NA removal")

    cpgs = 100
    expect_equal(dim(reduce_scMethrix(scm,n_cpg = cpgs,var="rand")),c(cpgs,ncol(scm)))
    expect_equal(dim(reduce_scMethrix(scm,n_cpg = cpgs,var="top")),c(cpgs,ncol(scm)))

    s1 <- reduce_scMethrix(scm,n_cpg = cpgs,var="top")
    s1 <- get_rowdata_stats(s1)
    scm <- get_rowdata_stats(scm)
    expect_equal(sort(rowData(s1)$sd,decreasing = TRUE),sort(rowData(scm)$sd,decreasing = TRUE)[1:cpgs])
    
  }))
})

test_that("dim_red_scMethrix", {

  expect_error(dim_red_scMethrix("not scMethrix"),msg.validateExp)
  expect_error(dim_red_scMethrix(scm.mem,type="not a type"),msg.validateArg)
  
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
