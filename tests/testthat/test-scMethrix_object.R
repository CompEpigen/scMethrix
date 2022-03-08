test_that("genome", {
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_equal(refGenome, unique(genome(scm)))
  }))
})

test_that("h5",{
  expect_false(S4Vectors::metadata(scm.mem)$is_h5) 
  expect_true(S4Vectors::metadata(scm.h5)$is_h5) 
})