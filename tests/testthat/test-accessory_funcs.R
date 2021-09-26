test_that("is_h5",{
  expect_error(is_h5("not a scMethrix"),"A valid scMethrix object needs to be supplied.")
  expect_true(is_h5(scm.h5))
  expect_false(is_h5(scm.mem))
})

test_that("has_cov",{
  expect_error(has_cov("not a scMethrix"),"A valid scMethrix object needs to be supplied.")
  expect_true(has_cov(scm.mem))
  expect_false(has_cov(remove_assay(scm.mem,assay="counts")))
})

test_that("get_sample_name", {
  expect_error(get_sample_name(5),"Must be a string file path")
  expect_equal("file.name",get_sample_name("c:/dir/dir.dir/file.name.extension"))
})

test_that("binarize", {
  expect_error(binarize("not numbers"),"non-numeric argument to binary operator")
  expect_equal(binarize(c(0,0,100,100,75,NA)),c(0,0,1,1,1,NA))
})

test_that("bin_granges",{
  expect_error(bin_granges("not granges"),"Input must be a Granges object")
  regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,100))
  expect_equal(length(bin_granges(regions,bin_size=10)),10) 
  expect_equal(reduce(bin_granges(regions,bin_size=10)),regions)
})

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

test_that("split_vector",{
  vec <- c(1,2,3,4,5,6,7,8)
  expect_error(split_vector(vec,size=1,chunks=1),"Invalid input. Must")
  expect_error(split_vector(vec,size="not numeric"),"Invalid input. Size")
  expect_error(split_vector(vec,percent="not numeric"),"Invalid input. Percent")
  expect_error(split_vector(vec,chunks="not numeric"),"Invalid input. Chunks")

  expect_equal(split_vector(vec,size=2),split_vector(vec,chunks=4))
  expect_equal(split_vector(vec,percent=25),split_vector(vec,chunks=4))
  expect_equal(unlist(split_vector(vec,percent=25)),vec)
  
  expect_equal(length(split_vector(vec,chunks=2.5)),3)
  expect_equal(length(split_vector(vec,percent=24)),5)
  expect_equal(length(split_vector(vec,size=2.5)),3)
  
  expect_equal(length(split_vector(vec,chunks=-10)),1)
  expect_equal(length(split_vector(vec,size=-10)),1)
  expect_equal(length(split_vector(vec,percent=-10)),1)
})

test_that("cast_granges",{
  expect_error(cast_granges("not a Granges"),"Invalid input class for regions. Must be a GRanges or data.frame-like")
  
  gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,100))
  expect_equal(gr,cast_granges(gr))
  
  df <- data.frame(chr=as.character(gr@seqnames),start=gr@ranges@start,end=gr@ranges@width)
  expect_equal(gr,cast_granges(df))
})

test_that("subset_ref_cpgs",{
  
  ref_cpgs = data.frame(chr="chr1",start=(1:5*2-1), end=(1:5*2))
  gen_cpgs = data.frame(chr="chr1",start=(6:10*2-1), end=(6:10*2))
  gen_cpgs = rbind(ref_cpgs[1:3,],gen_cpgs)
  rows <- sample(nrow(gen_cpgs))
  gen_cpgs <- gen_cpgs[rows, ] #randomize the order
  
  sub_cpgs = subset_ref_cpgs(ref_cpgs, gen_cpgs)
  
  expect_equal(ref_cpgs[1:3,],sub_cpgs)
  
})
