test_that("reduce_scMethrix", {
  
  expect_error(reduce_scMethrix("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  
    some_cpgs = nrow(scm)-1

    expect_error(reduce_scMethrix(scm="not scMethrix"),msg.validateExp)
    expect_error(reduce_scMethrix(scm,assay="not an assay"),msg.validateAssay)
    expect_error(reduce_scMethrix(scm,var="not a var"),msg.validateArg)
    expect_error(reduce_scMethrix(scm,n_cpg = "not an int"),msg.validateType)
    
    expect_true(identical.scm(scm,reduce_scMethrix(scm,n_cpg=nrow(scm))))
    
    expect_error(reduce_scMethrix(scm,n_cpg = 0,var="top"),"No CpGs left after reduction.")
    #s <- transform_assay(scm,new_assay="NA",trans = function(x) rep(NA,length(x)))
    #expect_error(reduce_scMethrix(s,n_cpg = some_cpgs, assay="NA",na.rm = T),"No CpGs left after reduction.")
    
    expect_equal(dim(reduce_scMethrix(scm,n_cpg = some_cpgs,var="rand")),c(some_cpgs,ncol(scm)))
    expect_equal(dim(reduce_scMethrix(scm,n_cpg = some_cpgs,var="top")),c(some_cpgs,ncol(scm)))

    s1 <- reduce_scMethrix(scm,n_cpg = some_cpgs,var="top")
    s1 <- get_rowdata_stats(s1)
    s2 <- get_rowdata_stats(scm)
    expect_equal(sort(rowData(s1)$sd,decreasing = TRUE),sort(rowData(s2)$sd,decreasing = TRUE)[1:some_cpgs])
    
  }))
})

test_that("dim_red_scMethrix", {

  expect_error(dim_red_scMethrix("not scMethrix"),msg.validateExp)
  expect_error(dim_red_scMethrix(scm.mem,type="not a type"),msg.validateArg)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {

    expect_error(dim_red_scMethrix(scm, assay = "score"),"Assay matrix cannot contain NAs")
    
    s <- suppressWarnings(impute_regions(scm))
    s <- reduce_scMethrix(s,assay="impute",n_cpg = 50)
    
    #PCA
    n_pc = 2
    pca <- dim_red_scMethrix(s, assay = "impute", n_pc = n_pc,type="PCA")
    expect_equal(reducedDimNames(pca),"PCA")
    expect_equal(dim(reducedDim(pca)),c(ncol(s),n_pc))
    expect_equal(length(pca@metadata$PCA_vars),n_pc)
    
    #UMAP
    umap <- dim_red_scMethrix(s, assay = "impute",type="UMAP")
    expect_equal(reducedDimNames(umap),"UMAP")
    expect_equal(dim(reducedDim(umap)),c(ncol(s),2))
    
    #tSNE
    tsne <- dim_red_scMethrix(s, assay = "impute",type="tSNE")
    expect_equal(reducedDimNames(tsne),"tSNE")
    expect_equal(dim(reducedDim(tsne)),c(ncol(s),2))
  }))
})
