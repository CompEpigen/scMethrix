
df1 <- data.table(chr=rep("chr1",5),start=1:5,end=2:6,value=0)
df2 <- data.table(chr=rep("chr1",5),start=3:7,end=4:8,value=0)
df3 <- data.table(chr=rep("chr1",5),start=6:10,end=7:11,value=0)

write.table(df1, file = file.path(tempdir(),'df1.bed'), row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df2, file = file.path(tempdir(),'df2.bed'), row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df3, file = file.path(tempdir(),'df3.bed'), row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)

files <- c("df1.bed","df2.bed","df3.bed")
files <- file.path(tempdir(),files)

test_that("read_index", {

  expect_equivalent(read_index(files),rbind(df1,df3)[,1:3])

})

test_that("read_bed_by_index", {
  
  index <- read_index(files)
  
  expect_equivalent(as.vector(read_bed_by_index(files[1],index)),c(rep(0,5),rep(NA,5)))
  expect_equivalent(as.vector(read_bed_by_index(files[3],index)),c(rep(NA,5),rep(0,5)))
  
})





