#---- genome -----------------------------------------------------------------------------------------------------------
test_that("genome", {
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_equal(refGenome, unique(genome(scm)))
  }))
})

#---- h5 ---------------------------------------------------------------------------------------------------------------
test_that("h5",{
  expect_false(S4Vectors::metadata(scm.mem)$is_h5) 
  expect_true(S4Vectors::metadata(scm.h5)$is_h5) 
})

#---- .validH5 ---------------------------------------------------------------------------------------------------------
test_that(".validH5",{
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_null(.validH5(scm))
    
    S4Vectors::metadata(scm)$is_h5 = "Not a bool"
    expect_false(is.null(.validH5(scm)))
    
    S4Vectors::metadata(scm)$is_h5 = FALSE
    expect_null(.validH5(scm))
    
    S4Vectors::metadata(scm)$is_h5 = NULL
    expect_false("is_h5" %in% names(S4Vectors::metadata(scm)))
    expect_false(is.null(.validH5(scm)))
  }))
})

#---- .validDims -------------------------------------------------------------------------------------------------------
test_that(".validDims",{
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    s <- scm
    expect_null(.validDims(s))
    s@assays@data[[1]] <- s@assays@data[[1]][,-1]
    expect_false(is.null(.validDims(s)))
    
    s <- scm
    expect_null(.validDims(s))
    s@assays@data[[1]] <- s@assays@data[[1]][-1,]
    expect_false(is.null(.validDims(s)))
  }))
})

#---- .validSamples ----------------------------------------------------------------------------------------------------
test_that(".validSamples",{
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_null(.validSamples(scm))
    newCols <- rep("Not a sample",ncol(scm))
    colnames(scm@assays@data[[1]]) <- newCols
    expect_null(.validSamples(scm))
    expect_identical(colnames(assays(scm,withDimnames = FALSE)[[1]]), newCols)
  }))
})

#---- .validFeatures ---------------------------------------------------------------------------------------------------
test_that(".validFeatures",{
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_null(.validFeatures(scm))
    newRows <- rep("Not a feature",nrow(scm))
    rownames(scm@assays@data[[1]]) <- newRows
    expect_null(.validFeatures(scm))
    expect_identical(rownames(assays(scm,withDimnames = FALSE)[[1]]), newRows)
  }))
})

#---- .validAssays -----------------------------------------------------------------------------------------------------
test_that(".validAssays",{
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_null(.validAssays(scm))
    scm@assays@data[[1]] <- as.data.frame(scm@assays@data[[1]])
    expect_false(is.null(.validAssays(scm)))
  }))
})

#---- .validRedDim -----------------------------------------------------------------------------------------------------
test_that(".validRedDim",{
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    #expect_null(.validRedDim(scm))
    # TODO: add this when function is fixed
  }))
})

#---- .validscMethrix --------------------------------------------------------------------------------------------------
test_that(".validscMethrix",{
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_null(.validscMethrix(scm))
  }))
})