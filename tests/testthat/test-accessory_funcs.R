#---- is_h5 ------------------------------------------------------------------------------------------------------------
test_that("is_h5",{
  expect_true(is_h5(scm.h5))
  expect_false(is_h5(scm.mem))
})

#---- has_cov ----------------------------------------------------------------------------------------------------------
test_that("has_cov",{
  expect_true(has_cov(scm.mem))
  expect_false(has_cov(remove_assay(scm.mem,assay="counts")))
})

#---- get_sample_name --------------------------------------------------------------------------------------------------
test_that("get_sample_name", {
  expect_error(get_sample_name(5),msg.validateType)
  expect_equal("file",get_sample_name("c:/dir/dir.dir/file"))
  expect_equal("file",get_sample_name("c:/dir/dir.dir/file.extension"))
  expect_equal("file",get_sample_name("c:/dir/dir.dir/file.extension.gz"))
  expect_equal("file.name",get_sample_name("c:/dir/dir.dir/file.name.extension.bz2"))
})

#---- normalize --------------------------------------------------------------------------------------------------------
test_that("normalize", {
  vals <- c(0,1,2,3,4,5)
  
  expect_equal(normalize(vals),vals/max(vals))
  expect_equal(normalize(vals, scale = TRUE),vals/max(vals))
  expect_equal(normalize(vals, min = min(vals), max = max(vals)*2),vals/max(vals)/2)
  expect_equal(normalize(vals, min = min(vals), max = max(vals), scale = TRUE),vals)
  expect_equal(normalize(vals, min = min(vals), max = max(vals)*2, scale = TRUE),vals*2)
})

#---- .pasteList -------------------------------------------------------------------------------------------------------
test_that(".pasteList", {

  vals = c("string1")
  expect_equal(.pasteList(vals),"'string1'")

  vals = c("string1","string2")
  expect_equal(.pasteList(vals),"'string1' and 'string2'")

  vals = c("string1","string2","string3")
  expect_equal(.pasteList(vals),"'string1', 'string2', and 'string3'")
})

#---- binarize ---------------------------------------------------------------------------------------------------------
test_that("binarize", {
  expect_error(binarize("not numbers"),"non-numeric argument to binary operator")
  expect_equal(binarize(c(0,0,100,100,75,NA)),c(0,0,1,1,1,NA))
})

#---- fill -------------------------------------------------------------------------------------------------------------
test_that("fill", {
  expect_equal(fill(c(0,0,100,100,75,NA)),c(0,0,100,100,75,0))
})

#---- start,split,stop_time --------------------------------------------------------------------------------------------
test_that("start,split,stop_time",{
 # expect_warning(stop_time())
#  expect_warning(split_time())
  
  start_time()
  
  Sys.sleep(1)
  elapsed_time <- capture.output(cat(split_time()))[1]
  elapsed_time <- as.numeric(strsplit(elapsed_time, "s")[[1]][1])
  expect_equal(elapsed_time, 1, tolerance = 0.1)
  
  Sys.sleep(2)
  
  elapsed_time <- capture.output(cat(split_time()))[1]
  elapsed_time <- as.numeric(strsplit(elapsed_time, "s")[[1]][1])
  expect_equal(elapsed_time, 2, tolerance = 0.25)
  
  Sys.sleep(3)
  
  elapsed_time <- capture.output(cat(stop_time()))[1]
  elapsed_time <- as.numeric(strsplit(elapsed_time, "s")[[1]][1])
  expect_equal(elapsed_time, 6, tolerance = 0.5)
  
})

#---- split_vector -----------------------------------------------------------------------------------------------------
test_that("split_vector",{
  vec <- c(1,2,3,4,5,6,7,8)
  expect_error(split_vector(vec,size=1,chunks=1),"Invalid input. Must contain 1 of")
  expect_error(split_vector(vec,size="not integer"),msg.validateType)
  expect_error(split_vector(vec,percent="not numeric"),msg.validateType)
  expect_error(split_vector(vec,chunks="not integer"),msg.validateType)

  expect_equal(split_vector(vec,size=2),split_vector(vec,chunks=4))
  expect_equal(split_vector(vec,percent=25),split_vector(vec,chunks=4))
  expect_equal(unlist(split_vector(vec,percent=25)),vec)
  
  expect_equal(length(split_vector(vec,chunks=3)),3)
  expect_equal(length(split_vector(vec,percent=24)),5)
  expect_equal(length(split_vector(vec,size=3)),3)
  
  expect_error(length(split_vector(vec,chunks=-10)),msg.validateValue)
  expect_error(length(split_vector(vec,size=-10)),msg.validateValue)
  expect_error(length(split_vector(vec,percent=-10)),msg.validateValue)
})


