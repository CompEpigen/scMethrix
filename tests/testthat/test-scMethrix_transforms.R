test_that("bin_scMethrix", {

  expect_error(bin_scMethrix("not scMethrix"))
  expect_error(bin_scMethrix(scm.h5))

  #invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  invisible(lapply(list(scm.mem), function(scm) {
    bin <- bin_scMethrix(scm, h5_dir = paste0(tempdir(),"/bin"), n_threads = 2)
    expect_equal(dim(bin),c(258,4))
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
      expect_equal(class(get_matrix(plus1,assay="plus1"))[[1]],"HDF5Matrix")
    }
    rm(plus1)
  }))
})

test_that("impute_by_iPCA", {
  imputation_test_helper(impute_by_iPCA)
})

test_that("impute_by_RF", {
  expect_warning(imputation_test_helper(impute_by_RF))
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
    expect_equal(nrow(set$training),57)
    expect_equal(nrow(set$test),229)
    
    expect_equivalent(scm,merge_scMethrix(set$training,set$test)) #TODO: validate this better
  }))
})

test_that("generate_random_subset", {
  
  expect_error(generate_training_set("not scMethrix"))
  
  #invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_warning(generate_random_subset(scm,n_cpgs = nrow(scm)+1))
    expect_warning(generate_random_subset(scm,n_cpgs = -1))
    
    m <- generate_random_subset(scm,n_cpgs = 100)
    expect_equal(nrow(m),100)
  }))
})


