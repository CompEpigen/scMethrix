df1 <- data.table(chr=rep("chr1",5),start=1:5,end=2:6,value=0)
df2 <- data.table(chr=rep("chr1",5),start=3:7,end=4:8,value=0)
df3 <- data.table(chr=rep("chr1",5),start=6:10,end=7:11,value=0)

files <- c("df1.bedgraph","df2.bedgraph","df3.bedgraph")
files <- file.path(tempdir(),files)

write.table(df1, file = files[1], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df2, file = files[2], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df3, file = files[3], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)

test_that("convert_HDF5_methrix", {

  expect_error(convert_HDF5_methrix("not scMethrix"))
  
  scm <- read_beds(files,h5=TRUE)

  expect_true(is_h5(scm))
  expect_equivalent(class(assays(scm)[[1]])[1],"DelayedMatrix")
  
  scm <- convert_HDF5_methrix(scm)
  
  expect_false(is_h5(scm))
  expect_equivalent(class(assays(scm)[[1]])[1],"matrix") 
  
})

test_that("convert_methrix", {
  
  expect_error(convert_methrix("not scMethrix"))
  
  scm <- read_beds(files,h5=FALSE)
  
  expect_false(is_h5(scm))
  expect_equivalent(class(assays(scm)[[1]])[1],"matrix") 
  
  scm <- convert_methrix(scm)
  
  expect_true(is_h5(scm))
  expect_equivalent(class(assays(scm)[[1]])[1],"HDF5Matrix")
  
})


