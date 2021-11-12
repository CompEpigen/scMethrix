test_that("is_h5",{
  expect_true(is_h5(scm.h5))
  expect_false(is_h5(scm.mem))
})

test_that("has_cov",{
  expect_true(has_cov(scm.mem))
  expect_false(has_cov(remove_assay(scm.mem,assay="counts")))
})

test_that("get_sample_name", {
  expect_error(get_sample_name(5),msg.validateType)
  expect_equal("file",get_sample_name("c:/dir/dir.dir/file"))
  expect_equal("file",get_sample_name("c:/dir/dir.dir/file.extension"))
  expect_equal("file",get_sample_name("c:/dir/dir.dir/file.extension.gz"))
  expect_equal("file.name",get_sample_name("c:/dir/dir.dir/file.name.extension.bz2"))
})

test_that("binarize", {
  expect_error(binarize("not numbers"),"non-numeric argument to binary operator")
  expect_equal(binarize(c(0,0,100,100,75,NA)),c(0,0,1,1,1,NA))
})

test_that("bin_granges",{
  expect_error(bin_granges(gr="not granges"),msg.validateType)
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

test_that("cast_granges",{
  expect_error(cast_granges("not a Granges"),"Invalid input class")
  
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

test_that(".validateExp",{
  expect_error(.validateExp("not scMethrix"),msg.validateExp)
  expect_true(.validateExp(scm.mem))
  expect_true(.validateExp(scm.h5))
})

test_that(".validateAssay",{
  expect_error(.validateAssay("not scMethrix"),msg.validateExp)
  expect_equivalent(.validateAssay(scm.mem,assay="score"),"score")
  expect_equivalent(.validateAssay(scm.mem,assay="sco"),"score")
  expect_error(.validateAssay(scm.mem,assay="not an assay"),msg.validateAssay)
})


test_that(".validateArg",{
  func <- function(var = c("banana","banjo")) {}
    
  var = "banana"
  expect_equivalent(.validateArg(var,func),"banana")    
  var = "banjo"
  expect_equivalent(.validateArg(var,func),"banjo")  
  var = "bAnA"
  expect_equivalent(.validateArg(var,func),"banana")  
  expect_error(.validateArg(var,func,ignore.case = F), msg.validateArg)  
  var = c("banana","banjo")
  expect_equivalent(.validateArg(var,func),"banana")    
  var = "ban"
  expect_error(.validateArg(var,func), msg.validateArg)  
  var = "bad input"
  expect_error(.validateArg(var,func), msg.validateArg) 
    
  #TODO: test input for argument list
})

test_that(".validateType",{
  
  expect_error(.validateType(input = "an input"),           "No valid type specified")
  expect_error(.validateType(input = "an input",            type = "not a type"),msg.validateArg)
  
  expect_true (.validateType(input = 10,                     type = "integer"))
 # expect_true (.validateType(input = c(10,20,30),            type = "integer"))
  expect_true (.validateType(input = 10,                     type = "INT"))
  expect_error(.validateType(input = "not an int",           type = "integer"),msg.validateType)
  #expect_error(.validateType(input = list(10,"not an int"),  type = "integer"),msg.validateType)
  
  expect_true (.validateType(input = 10.5,                   type = "numeric"))
  #expect_true (.validateType(input = c(10.5,20.5,30.5),      type = "numeric"))
  expect_error(.validateType(input = "not an num",           type = "numeric"),msg.validateType)
  #expect_error(.validateType(input = list(10.5,"not an num"),type = "numeric"),msg.validateType)
  
  expect_true (.validateType(input = "A",                    type = "character"))
  #expect_true (.validateType(input = c("A","B"),             type = "character"))
  expect_error(.validateType(input = "not a char",           type = "character"),msg.validateType)
  #expect_error(.validateType(input = c("A","not a char"),    type = "character"),msg.validateType)
  
  expect_true (.validateType(input = "str1",                 type = "string"))
  #expect_true (.validateType(input = c("str1","str2"),       type = "string"))
  expect_error(.validateType(input = 0,                      type = "string"),msg.validateType)
  #expect_error(.validateType(input = list("str1",0),         type = "string"),msg.validateType)
  
  expect_true (.validateType(input = c(1,2,3),               type = "vector"))
  #expect_error(.validateType(input = 1,                      type = "vector"),msg.validateType)
  
  expect_true (.validateType(input = list(1,2),              type = "list"))
 # expect_true (.validateType(input = c(list(1,2),list(1,2)), type = "list"))
  expect_error(.validateType(input = "not lst",              type = "list"),msg.validateType)
  #expect_true (.validateType(input = c(list(1,2),"not lst"), type = "list"),msg.validateType)
  
  expect_true (.validateType(input = TRUE,                   type = "boolean"))
  #expect_true (.validateType(input = c(TRUE,TRUE),           type = "boolean"))
  expect_error(.validateType(input = "not a bool",           type = "boolean"),msg.validateType)  
  #expect_error(.validateType(input = c(T,"not a bool"),      type = "boolean"),msg.validateType) 
  expect_true (.validateType(input = TRUE,                   type = "logical"))
  expect_error(.validateType(input = "not a bool",           type = "logical"),msg.validateType)
  
  tmpfile <- tempfile()
  file.create(tmpfile)
  expect_true (.validateType(input = tmpfile,                type = "file"))
  #expect_true (.validateType(input = c(tmpfile,tmpfile),     type = "file"))
  expect_error(.validateType(input = "not a file",           type = "file"),msg.validateType)
  
  expect_true (.validateType(input = tempdir(),              type = "directory"))
  expect_error(.validateType(input = "not an directory",     type = "directory"),msg.validateType)
  
  gr <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), IRanges(1:10, width=10:1))
  expect_true (.validateType(input = GRanges(),              type = "GRanges"))
  expect_true (.validateType(input = gr,                     type = "GRanges"))
  #expect_true (.validateType(input = c(gr,gr),               type = "GRanges"))
  expect_error(.validateType(input = "not an Granges",       type = "GRanges"),msg.validateType)
  
  expect_true (.validateType(input = sum,                    type = "function"))
  #expect_true (.validateType(input = c(sum,sum),             type = "function"))
  #expect_true (.validateType(input = function(x) x+1,        type = "function"))
  expect_error(.validateType(input = "not a function",       type = "function"),msg.validateType)
  
  expect_true (.validateType(input = NULL,                   type = "null"))
  expect_error(.validateType(input = "not null",             type = "null"),msg.validateType)
  
})
