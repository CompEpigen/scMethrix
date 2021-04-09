
df1 <- data.table(chr=rep("chr1",5),start=1:5,end=2:6,value=0)
df2 <- data.table(chr=rep("chr1",5),start=3:7,end=4:8,value=0)
df3 <- data.table(chr=rep("chr1",5),start=6:10,end=7:11,value=0)
df4 <- data.table(chr=rep("chr2",5),start=1:10,end=2:11,value=0)

files <- c("df1.bedgraph","df2.bedgraph","df3.bedgraph","df4.bedgraph")
files <- file.path(tempdir(),files)

write.table(df1, file = files[1], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df2, file = files[2], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df3, file = files[3], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df4, file = files[4], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)

scm.h5 <- read_beds(files,h5=TRUE)
scm.mem <- read_beds(files,h5=FALSE)


test_that("convert_HDF5_methrix", {

  expect_error(convert_HDF5_methrix("not scMethrix"))
  
  expect_true(is_h5(scm.h5))
  expect_equivalent(class(get_matrix(scm.h5))[1],"DelayedMatrix")
  
  scm <- convert_HDF5_methrix(scm.h5)
  
  expect_false(is_h5(scm))
  expect_equivalent(class(get_matrix(scm))[1],"matrix") 
  rm(scm)
  
})

test_that("convert_methrix", {
  
  expect_error(convert_methrix("not scMethrix"))
  
  expect_false(is_h5(scm.mem))
  expect_equivalent(class(get_matrix(scm.mem))[1],"matrix") 
  
  scm <- convert_methrix(scm.mem)
  
  expect_true(is_h5(scm))
  expect_equivalent(class(get_matrix(scm))[1],"HDF5Matrix")
  rm(scm)
  
})

test_that("subset_scMethrix", {
  
  expect_error(subset_scMethrix("not scMethrix"))
  expect_warning(subset_scMethrix(scm.h5))
  
  samples <- c("df1","df3")
  s <- subset_scMethrix(scm.h5, samples = samples)
  expect_equivalent(dim(s),c(20,2))
  
  contigs <- c("chr1")
  s <- subset_scMethrix(scm.h5, contigs = contigs)
  expect_equivalent(dim(s),c(10,4))
  
  regions <- GRanges(seqnames = "chr1", ranges = IRanges(1,5)) 
  s <- subset_scMethrix(scm.h5, regions = regions)
  expect_equivalent(dim(s),c(5,4))
  
  s <- subset_scMethrix(scm.h5, samples = samples, contigs = contigs, regions = regions)
  expect_equivalent(dim(s),c(5,2))
  
})

test_that("region_filter", {
  
  expect_error(region_filter("not scMethrix"))
  expect_warning(region_filter(scm.h5))
  
  regions <- GRanges(seqnames = "chr1", ranges = IRanges(1,5)) 
  s <- region_filter(scm.h5, regions = regions)
  expect_equivalent(dim(s),c(16,4))
  
})

test_that("get_matrix", {

  expect_error(get_matrix("not scMethrix"))
  expect_warning(get_matrix(scm.h5,add_loci=FALSE, in_granges = TRUE))

  m <- get_matrix(scm.h5)
  expect_equivalent(dim(m),c(20,4))  
  expect_equivalent(class(m)[1],"DelayedMatrix")
  
  m <- get_matrix(scm.h5,add_loci=TRUE)
  expect_equivalent(dim(m),c(20,7))  
  expect_equivalent(class(m)[1],"data.table")
  
  m <- get_matrix(scm.h5,add_loci=TRUE, in_granges = TRUE)
  expect_equivalent(seqnames(m)@lengths,c(10,10))
  expect_equivalent(dim(mcols(m)),c(20,4))
  expect_equivalent(class(m)[1],"GRanges")

  m <- get_matrix(scm.mem)
  expect_equivalent(dim(m),c(20,4))  
  expect_equivalent(class(m)[1],"matrix")
  
  m <- get_matrix(scm.mem,add_loci=TRUE)
  expect_equivalent(dim(m),c(20,7))  
  expect_equivalent(class(m)[1],"data.table")
  
  m <- get_matrix(scm.mem,add_loci=TRUE, in_granges = TRUE)
  expect_equivalent(seqnames(m)@lengths,c(10,10))
  expect_equivalent(dim(mcols(m)),c(20,4))
  expect_equivalent(class(m)[1],"GRanges")
  
})


test_that("remove_uncovered", {
  
  expect_error(get_matrix("not scMethrix"))

  h5 <- remove_uncovered(subset_scMethrix(scm.h5,samples="df1"))
  expect_equivalent(dim(h5),c(5,1))

  mem <- remove_uncovered(subset_scMethrix(scm.mem,samples="df1"))
  expect_equivalent(dim(mem),c(5,1))
  
  expect_equivalent(h5,mem)
  
  rm(h5,mem)
  
})

test_that("get_stats", {
  
  expect_error(get_matrix("not scMethrix"))
  
  h5 = get_stats(scm.h5)
  mem = get_stats(scm.mem)

  expect_equivalent(dim(h1),c(8,5))
  expect_equivalent(dim(h2),c(4,4))
  expect_equivalent(h1,get_stats(scm.h5,per_chr=FALSE))
  expect_equivalent(h2,get_stats(scm.mem,per_chr=FALSE))
  
  
  rm(h5,mem)
})


