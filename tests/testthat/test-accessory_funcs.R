
df1 <- data.table(chr=rep("chr1",5),start=1:5,end=2:6,value=0)
df2 <- data.table(chr=rep("chr1",5),start=3:7,end=4:8,value=0)
df3 <- data.table(chr=rep("chr1",5),start=6:10,end=7:11,value=0)

files <- c("df1.bedgraph","df2.bedgraph","df3.bedgraph")
files <- file.path(tempdir(),files)

write.table(df1, file = files[1], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df2, file = files[2], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df3, file = files[3], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)

scm.h5 <- read_beds(files,h5=TRUE)
scm.mem <- read_beds(files,h5=FALSE)

test_that("is_ondisk",{
  expect_true(is_h5(scm.h5))
  expect_false(is_h5(scm.mem))
})

test_that("get_sample_name", {
  expect_equal("file.name",get_sample_name("c:/dir/dir.dir/file.name.extension"))
})

test_that("chunk_granges",{
  
  expect_error(chunk_granges("not scMethrix"))
  expect_error(chunk_granges("not scMethrix"))
  
  expect_true(all(rowRanges(scm.h5) == unlist(chunk_granges(rowRanges(scm.h5),factor=factor,percent=percent,num=num))))
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





test_that("cast_granges",{})

test_that("order_by_sd",{})

test_that("subset_scMethrix",{})

test_that("get_stats",{})