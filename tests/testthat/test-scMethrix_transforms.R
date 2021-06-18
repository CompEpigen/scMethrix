test_that("bin_scMethrix", {
  
  expect_error(bin_scMethrix("not scMethrix"))
  expect_error(bin_scMethrix(scm.h5))
  
  #invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  invisible(lapply(list(scm.mem), function(scm) {
    bin <- bin_scMethrix(scm,bin_size = 10000000, h5_dir = paste0(tempdir(),"/bin"))
    expect_equivalent(dim(bin),c(77,4))
    rm(bin)
  }))
})

test_that("transform_assay", {
  
  expect_error(transform_assay("not scMethrix"))
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(transform_assay(scm,trans="not closure"))
    expect_warning(transform_assay(scm,trans=function(x) x+1,assay="score",new_assay="score"))
    plus1 <- transform_assay(scm,trans=function(x) x+1,assay="score",new_assay="plus1")
    expect_false(isTRUE(all.equal(assays(scm), assays(plus1))))
    expect_equivalent(get_matrix(scm)+1,get_matrix(plus1,assay="plus1"))
    if (is_h5(scm)) {
      expect_equivalent(class(get_matrix(plus1,assay="plus1"))[[1]],"HDF5Matrix")
    }
    rm(plus1)
  }))
})

test_that("impute_by_iPCA", {
  imputation_test_helper(impute_by_iPCA)
})

test_that("impute_by_RF", {
  imputation_test_helper(impute_by_RF)
})

test_that("impute_by_kNN", {
  imputation_test_helper(impute_by_kNN)
})

test_that("generate_training_set", {
  
  expect_error(generate_training_set("not scMethrix"))
  
  #invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(generate_training_set(scm,training_prop = 2))
    
    set <- generate_training_set(scm,training_prop = 0.2)
    expect_equivalent(nrow(set$training),20)
    expect_equivalent(nrow(set$test),80)
    
    expect_equivalent(scm,merge_scMethrix(set$training,set$test)) #TODO: validate this better
  }))
})



