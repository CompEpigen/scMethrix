test_that("bin_scMethrix", {

  expect_error(bin_scMethrix("not scMethrix"),"A valid scMethrix object needs to be supplied")
  expect_error(bin_scMethrix(scm.h5),"Output directory must be specified")

  #invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  invisible(lapply(list(scm.mem), function(scm) {
    bin <- bin_scMethrix(scm, h5_dir = paste0(tempdir(),"/bin"), n_threads = 2)
    expect_equal(dim(bin),c(258,4))
    rm(bin)
  }))
})

test_that("transform_assay", {
  
  expect_error(transform_assay("not scMethrix"),"A valid scMethrix object needs to be supplied")
  trans <- function(x) x+1
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(transform_assay(scm,trans="not closure"),"A valid transform function must be specified")
    expect_warning(transform_assay(scm, trans=trans, assay="score",new_assay="score"))
    plus1 <- transform_assay(scm, trans=trans, assay="score",new_assay="plus1")
    expect_false(isTRUE(all.equal(assays(scm), assays(plus1))))
    expect_equivalent(get_matrix(scm)+1,get_matrix(plus1,assay="plus1"))
    if (is_h5(scm)) {
      expect_equal(class(get_matrix(plus1,assay="plus1"))[[1]],"HDF5Matrix")
    }
    rm(plus1)
  }))
})

test_that("impute_regions", {
  expect_error(impute_regions("not scMethrix"),"A valid scMethrix object needs to be supplied")
  
  suppressWarnings(
    lapply(list("kNN","iPCA","RF"), function (method) {
    #invisible(lapply(list(scm.mem,scm.h5), function(scm) {
      invisible(lapply(list(scm.mem), function(scm) {
        expect_error(impute_regions(scm,assay = "not an assay"))
        expect_error(impute_regions(scm,new_assay = "score"))
        expect_warning(impute_regions(scm,new_assay = "counts",type=method))
        
        impute = impute_regions(scm,new_assay="impute",type=method)
        expect_true("impute" %in% SummarizedExperiment::assayNames(impute))
        
        sco <- get_matrix(impute,assay="score")
        imp <- get_matrix(impute,assay="impute")
        NAs <- which(is.na(sco))
        nonNAs <- which(!is.na(sco))
        
        expect_true(anyNA(sco) && !anyNA(imp))
        expect_equivalent(sco[nonNAs],imp[nonNAs])
        expect_false(all(sco[NAs] %in% imp[NAs]))
      }))
    })
  )
})

test_that("generate_training_set", {
  
  expect_error(generate_training_set("not scMethrix"),"A valid scMethrix object needs to be supplied")
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(generate_training_set(scm,training_prop = 2),"training_prop must in the range of")
    
    set <- generate_training_set(scm,training_prop = 0.2)
    expect_equal(nrow(set$training),57)
    expect_equal(nrow(set$test),229)
    
    expect_equivalent(scm,merge_scMethrix(set$training,set$test)) #TODO: validate this better
  }))
})

test_that("generate_random_subset", {
  
  expect_error(generate_training_set("not scMethrix"),"A valid scMethrix object needs to be supplied")
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_warning(generate_random_subset(scm,n_cpgs = nrow(scm)+1))
    expect_warning(generate_random_subset(scm,n_cpgs = -1))
    
    m <- generate_random_subset(scm,n_cpgs = 100)
    expect_equal(nrow(m),100)
  }))
})