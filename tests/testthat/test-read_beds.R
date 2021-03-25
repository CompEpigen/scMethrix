
df1 <- data.table(chr=rep("chr1",5),start=1:5,end=2:6,value=0)
df2 <- data.table(chr=rep("chr1",5),start=3:7,end=4:8,value=0)
df3 <- data.table(chr=rep("chr1",5),start=6:10,end=7:11,value=0)

files <- c("df1.bedgraph","df2.bedgraph","df3.bedgraph")
files <- file.path(tempdir(),files)

write.table(df1, file = files[1], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df2, file = files[2], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df3, file = files[3], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)

test_that("read_index", {

  expect_equivalent(read_index(files),rbind(df1,df3)[,1:3])

})

test_that("read_bed_by_index", {
  
  index <- read_index(files)
  
  expect_equivalent(as.vector(read_bed_by_index(files[1],index)),c(rep(0,5),rep(NA,5)))
  expect_equivalent(as.vector(read_bed_by_index(files[3],index)),c(rep(NA,5),rep(0,5)))
  
})

test_that("read_bed - Input errors", {
  
  expect_error(read_beds(NULL))
  expect_error(read_beds(tempfile))
  
})

test_that("read_bed - HDF5", {
  
  path <- file.path(tempdir(),paste0("sse-",sample.int(10000, 1)))
  
  scm1 <- read_beds(files,h5=TRUE,h5_dir=path)
  scm2 <- load_HDF5_scMethrix(dir=path)
    
  expect_equivalent(class(scm1)[1],class(scm2)[1],"scMethrix")
  expect_equivalent(class(get_matrix(scm1))[[1]],class(get_matrix(scm2))[[1]],"DelayedMatrix")
  expect_equivalent(dim(scm1),dim(scm2),c(10,3))
  
})

test_that("read_bed - in-memory", {
  
  scm <- read_beds(files,h5=FALSE)
  
  expect_equivalent(class(scm)[1],"scMethrix")
  expect_equivalent(class(get_matrix(scm))[[1]],"matrix")
  expect_equivalent(dim(scm),c(10,3))
  
})

test_that("read_bed - HDF5 and in-memory equivalence", {
  
  scm.hdf <- read_beds(files,h5=TRUE)
  scm.mem <- read_beds(files,h5=FALSE)
  
  expect_equivalent(as.matrix(assays(scm.hdf)$score),assays(scm.mem)$score)
  expect_equivalent(rowRanges(scm.hdf),rowRanges(scm.mem))
})


