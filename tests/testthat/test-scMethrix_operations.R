test_that("convert_HDF5_methrix", {

  scm <- NULL #load the h5 object

  expect_true(is_h5(scm))
  expect_equivalent(class(assays(scm)[[1]])[1],"DelayedMatrix")
  
  scm <- convert_HDF5_methrix(scm)
  
  expect_false(is_h5(scm))
  expect_equivalent(typeof(assays(scm)[[1]])[1],"matrix")
  
})

test_that("convert_methrix", {
  
  scm <- NULL #load the h5 object
  
  expect_false(is_h5(scm))
  expect_equivalent(typeof(assays(scm)[[1]])[1],"matrix")
  
  scm <- convert_methrix(scm)
  
  expect_true(is_h5(scm))
  expect_equivalent(class(assays(scm)[[1]])[1],"DelayedMatrix")
  
})


