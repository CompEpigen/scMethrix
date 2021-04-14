test_that("is_ondisk",{
  expect_true(is_h5(scm.h5))
  expect_false(is_h5(scm.mem))
})

test_that("get_sample_name", {
  expect_equal("file.name",get_sample_name("c:/dir/dir.dir/file.name.extension"))
})

test_that("split_granges",{
  
  expect_error(split_granges("not scMethrix"))
  expect_error(split_granges(rowRanges(scm.h5)))
  expect_error(split_granges(rowRanges(scm.h5),factor=2,num=3))
  
  factor = 2
  percent = 50 
  num = 2
  
  expect_true(all(rowRanges(scm.h5) == unlist(split_granges(rowRanges(scm.h5),factor=factor))))
  expect_true(all(rowRanges(scm.h5) == unlist(split_granges(rowRanges(scm.h5),percent=percent))))
  expect_true(all(rowRanges(scm.h5) == unlist(split_granges(rowRanges(scm.h5),num=num))))
})

test_that("start,split,stop_time",{
  
  expect_error(split_time())
  expect_error(stop_time())
  
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
  
  expect_error(split_time())
  expect_error(stop_time())
  
})

test_that("split_vector",{
  
  vec <- c(1,2,3,4,5,6,7,8)
  expect_equivalent(split_vector(vec,4,by="size"),split_vector(vec,2,by="chunk"))
  
})

test_that("cast_granges",{})

test_that("order_by_sd",{})

test_that("get_stats",{})