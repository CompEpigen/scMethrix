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

test_that("convert_HDF5_methrix", {

  expect_error(convert_HDF5_methrix("not scMethrix"))
  
  expect_true(is_h5(scm.h5))
  expect_equivalent(class(assays(scm.h5)[[1]])[1],"DelayedMatrix")
  
  scm <- convert_HDF5_methrix(scm.h5)
  
  expect_false(is_h5(scm))
  expect_equivalent(class(assays(scm)[[1]])[1],"matrix") 
  rm(scm)
  
})

test_that("convert_methrix", {
  
  expect_error(convert_methrix("not scMethrix"))
  
  expect_false(is_h5(scm.mem))
  expect_equivalent(class(assays(scm.mem)[[1]])[1],"matrix") 
  
  scm <- convert_methrix(scm.mem)
  
  expect_true(is_h5(scm))
  expect_equivalent(class(assays(scm)[[1]])[1],"HDF5Matrix")
  rm(scm)
  
})



test_that("subset_scMethrix", {
  
  
  

  
})

test_that("get_matrix") {

  expect_error(get_matrix("not scMethrix"))

  m <- get_matrix(scm.h5)
  expect_equivalent(dim(m),c(10,3))  
  expect_equivalent(class(m)[1],"data.frame")
  
  m <- get_matrix(scm.h5,add_loci=TRUE)
  expect_equivalent(dim(m),c(10,6))  
  expect_equivalent(class(m)[1],"data.frame")
  
  m <- get_matrix(scm.h5,add_loci=TRUE, in_granges = TRUE)
  expect_equivalent(seqnames(m)@lengths,10)
  expect_equivalent(dim(mcols(m)),c(10,3))
  expect_equivalent(class(m)[1],"GRanges")

  expect_warning(get_matrix(scm.h5,add_loci=FALSE, in_granges = TRUE))
  
}



