test_that("is_ondisk",{
  expect_true(is_h5(scm.h5))
  expect_false(is_h5(scm.mem))
})

test_that("get_sample_name", {
  expect_equal("file.name",get_sample_name("c:/dir/dir.dir/file.name.extension"))
})

test_that("binarize", {
  expect_equal(binarize(75),1)
  expect_equal(binarize(25),0)
})

test_that("bin_granges",{
  regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,100))
  expect_equivalent(length(bin_granges(regions,bin_size=10)),10) 
  expect_equivalent(reduce(bin_granges(regions,bin_size=10)),regions)
})

test_that("split_granges",{
  regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,100))
  regions <- bin_granges(regions,bin_size=1)
  
  expect_error(split_granges("not granges"))
  expect_error(split_granges(regions))
  expect_error(split_granges(regions,chunks=2,size=3))
  
  chunks = 3
  percent = 33
  size = 33
  
  expect_equivalent(length(split_granges(regions,chunks=chunks)),3)
  expect_equivalent(length(split_granges(regions,percent=percent)),4)
  expect_equivalent(length(split_granges(regions,percent=percent)[[4]]),1)
  expect_equivalent(split_granges(regions,percent=percent),split_granges(regions,size=size))
  expect_true(all(regions == unlist(split_granges(regions,size=size))))

})

test_that("start,split,stop_time",{
  expect_warning(split_time())
  expect_warning(stop_time())
  
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

test_that("split_vector",{
  vec <- c(1,2,3,4,5,6,7,8)
  expect_equivalent(split_vector(vec,4,by="size"),split_vector(vec,2,by="chunk"))
})

test_that("cast_granges",{
  expect_error(cast_granges("not a Granges"))
  
  gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,100))
  expect_equivalent(gr,cast_granges(gr))
  
  df <- data.frame(chr=as.character(gr@seqnames),start=gr@ranges@start,end=gr@ranges@width)
  expect_equivalent(gr,cast_granges(df))
})

test_that("subset_ref_cpgs",{
  
  ref_cpgs = data.frame(chr="chr1",start=(1:5*2-1), end=(1:5*2))
  gen_cpgs = data.frame(chr="chr1",start=(6:10*2-1), end=(6:10*2))
  gen_cpgs = rbind(ref_cpgs[1:3,],gen_cpgs)
  rows <- sample(nrow(gen_cpgs))
  gen_cpgs <- gen_cpgs[rows, ] #randomize the order
  
  sub_cpgs = subset_ref_cpgs(ref_cpgs, gen_cpgs)
  
  expect_equivalent(ref_cpgs[1:3,],sub_cpgs)
  
})
