data('scMethrix_data')
mem <- NULL
tbx <- NULL

test_that("is_ondisk",{
  expect_true(is_ondisk(tbx))
  expect_false(is_ondisk(mem))
})

test_that("get_sample_name", {
  expect_equal("filename",get_sample_name("c:/dir/dir.dir/filename.extension.extension"))
})

test_that("get_files",{
  
})

test_that("chunk_granges",{
  expect_true(all(rowRanges(mem) == unlist(chunk_granges(rowRanges(mem),factor=factor,percent=percent,num=num))))
  expect_true(all(rowRanges(tbx) == unlist(chunk_granges(rowRanges(tbx),factor=factor,percent=percent,num=num))))
})

test_that("cast_granges",{})

test_that("order_by_sd",{})

test_that("subset_scMethrix",{})

test_that("get_stats",{})