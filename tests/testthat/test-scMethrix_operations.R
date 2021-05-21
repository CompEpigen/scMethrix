
test_that("transform_assay", {
  invisible(lapply(list(scm.mem,scm.h5), function(scm) { 
    expect_error(transform_assay(m="not scMethrix"))
    expect_error(transform_assay(scm,trans="not closure"))
    expect_error(transform_assay(scm,trans=function(x) x+1,assay="not an assay"))
    expect_error(transform_assay(scm,trans=function(x) x+1,assay="score",name="score"))
    plus1 <- transform_assay(scm,trans=function(x) x+1,assay="score",name="plus1")
    expect_false(isTRUE(all.equal(assays(scm), assays(plus1))))
    expect_equivalent(get_matrix(scm)+1,get_matrix(plus1,type="plus1"))
    rm(plus1)
  }))
})

test_that("remove_assay", {
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(remove_assay(m="not scMethrix"))
    expect_error(remove_assay(scm, assay="not an assay"))
    expect_error(remove_assay(scm, assay="score"))
    plus1 <- transform_assay(scm,trans=function(x) x+1,assay="score",name="plus1")
    expect_false(isTRUE(all.equal(assays(scm), assays(plus1))))
    plus1 <- remove_assay(plus1, assay="plus1")
    expect_true(all.equal(assays(scm), assays(plus1)))
  }))
})

test_that("merge_scMethrix", {
    expect_error(merge_scMethrix(m1=scm.mem,m2="not scMethrix"))
  
    expect_error(merge_scMethrix(scm.mem,scm.mem,by="col")) #same samples
    expect_error(merge_scMethrix(scm.mem[1,1],scm.mem[2,2],by="col")) #different regions
    
    expect_error(merge_scMethrix(scm.mem[,1],scm.mem[,2],by="row")) #different samples
    expect_error(merge_scMethrix(scm.mem,scm.mem,by="row")) #same regions
    
    s <- scm.mem
    assays(s) <- assays(s)["score"]
    expect_warning(merge_scMethrix(s[,1],scm.mem[,2],by="col")) #different assays
    
    expect_equivalent(merge_scMethrix(scm.mem[1],scm.mem[2:nrow(scm.mem)],by="row"),scm.mem)
    expect_equivalent(merge_scMethrix(scm.mem[,1],scm.mem[,2:ncol(scm.mem)],by="col"),scm.mem)
})

test_that("convert_HDF5_scMethrix", {

  expect_error(convert_HDF5_scMethrix("not scMethrix"))
  
  expect_true(is_h5(scm.h5))
  expect_equivalent(class(get_matrix(scm.h5))[1],"HDF5Matrix")
  
  scm <- convert_HDF5_scMethrix(scm.h5)
  
  expect_false(is_h5(scm))
  expect_equivalent(class(get_matrix(scm))[1],"matrix") 
  rm(scm)
  
})

test_that("convert_scMethrix", {
  
  expect_error(convert_scMethrix("not scMethrix"))
  
  expect_false(is_h5(scm.mem))
  expect_equivalent(class(get_matrix(scm.mem))[1],"matrix") 
  
  scm <- convert_scMethrix(scm.mem)
  
  expect_true(is_h5(scm))
  expect_equivalent(class(get_matrix(scm))[1],"HDF5Matrix")
  rm(scm)
  
})

test_that("subset_scMethrix", {
  
  expect_error(subset_scMethrix("not scMethrix"))
  expect_warning(subset_scMethrix(scm.h5))
  
  samples <- c("C1","C3")
  contigs <- c("chr1","chr2")
  regions <- GRanges(seqnames = c("chr1","chr2","chr3","chr4"), ranges = IRanges(1,100000000)) 
  
  # Subset H5 by include
  
  s <- subset_scMethrix(scm.h5, samples = samples, by="include")
  expect_equivalent(dim(s),c(100,2))
  expect_equivalent(samples,colData(s)@rownames)

  s <- subset_scMethrix(scm.h5, contigs = contigs, by="include")
  expect_equivalent(dim(s),c(9,4))
  expect_equivalent(contigs,as.character(seqnames(s)@values))
  
  s <- subset_scMethrix(scm.h5, regions = regions, by="include")
  expect_equivalent(dim(s),c(10,4))

  s <- subset_scMethrix(scm.h5, samples = samples, contigs = contigs, regions = regions, by="include")
  expect_equivalent(dim(s),c(4,2))

  # Subset mem by include
  
  s <- subset_scMethrix(scm.mem, samples = samples, by="include")
  expect_equivalent(dim(s),c(100,2))

  s <- subset_scMethrix(scm.mem, contigs = contigs, by="include")
  expect_equivalent(dim(s),c(9,4))
  
  s <- subset_scMethrix(scm.mem, regions = regions, by="include")
  expect_equivalent(dim(s),c(10,4))
  
  s <- subset_scMethrix(scm.mem, samples = samples, contigs = contigs, regions = regions, by="include")
  expect_equivalent(dim(s),c(4,2))
  
  # Subset H5 by exclude
  
  s <- subset_scMethrix(scm.h5, samples = samples, by = "exclude")
  expect_equivalent(dim(s),c(100,2))
  expect_equivalent(length(intersect(colData(s)$colData,samples)),0)
  
  s <- subset_scMethrix(scm.h5, contigs = contigs, by = "exclude")
  expect_equivalent(dim(s),c(91,4))
  #expect_equivalent(contigs,as.character(seqnames(s)@values))
  
  s <- subset_scMethrix(scm.h5, regions = regions, by = "exclude")
  expect_equivalent(dim(s),c(90,4))
  
  s <- subset_scMethrix(scm.h5, samples = samples, contigs = contigs, regions = regions, by = "exclude")
  expect_equivalent(dim(s),c(85,2))
  
  # Subset mem by exclude
  
  s <- subset_scMethrix(scm.mem, samples = samples, by = "exclude")
  expect_equivalent(dim(s),c(100,2))
  
  s <- subset_scMethrix(scm.mem, contigs = contigs, by = "exclude")
  expect_equivalent(dim(s),c(91,4))
  
  s <- subset_scMethrix(scm.mem, regions = regions, by = "exclude")
  expect_equivalent(dim(s),c(90,4))
  
  s <- subset_scMethrix(scm.mem, samples = samples, contigs = contigs, regions = regions, by = "exclude")
  expect_equivalent(dim(s),c(85,2))
  
})

test_that("get_matrix", {

  expect_error(get_matrix("not scMethrix"))
  expect_warning(get_matrix(scm.h5,add_loci=FALSE, in_granges = TRUE))

  m <- get_matrix(scm.h5)
  expect_equivalent(dim(m),c(100,4))  
  expect_equivalent(class(m)[1],"HDF5Matrix")
  
  m <- get_matrix(m=scm.h5,add_loci=TRUE)
  expect_equivalent(dim(m),c(100,7))  
  expect_equivalent(class(m)[1],"data.table")
  
  m <- get_matrix(scm.h5,add_loci=TRUE, in_granges = TRUE)
 # expect_equivalent(seqnames(m)@lengths,c(10,8))
  expect_equivalent(dim(mcols(m)),c(100,4))
  expect_equivalent(class(m)[1],"GRanges")

  m <- get_matrix(scm.mem)
  expect_equivalent(dim(m),c(100,4))  
  expect_equivalent(class(m)[1],"matrix")
  
  m <- get_matrix(scm.mem,add_loci=TRUE)
  expect_equivalent(dim(m),c(100,7))  
  expect_equivalent(class(m)[1],"data.table")
  
  m <- get_matrix(scm.mem,add_loci=TRUE, in_granges = TRUE)
  #expect_equivalent(seqnames(m)@lengths,c(10,8))
  expect_equivalent(dim(mcols(m)),c(100,4))
  expect_equivalent(class(m)[1],"GRanges")
  
})


test_that("remove_uncovered", {
  
  expect_error(get_matrix("not scMethrix"))

  h5 <- subset_scMethrix(scm.h5,samples=c("C1","C2"),by = "include")
  expect_equivalent(dim(h5),c(100,2))
  expect_equivalent(dim(remove_uncovered(h5)),c(81,2))

  mem <- subset_scMethrix(scm.mem,samples=c("C1","C2"),by = "include")
  expect_equivalent(dim(h5),c(100,2))
  expect_equivalent(dim(remove_uncovered(mem)),c(81,2))
  
  expect_equivalent(h5,mem)
  
  rm(h5,mem)
  
})

test_that("get_stats", {
  
  expect_error(get_matrix("not scMethrix"))

  
  chr <- length(seqlengths(rowRanges(scm.h5)))
  samples <- nrow(colData(scm.h5))
  
  expect_equivalent(dim(get_stats(scm.h5)),c(chr*samples,5))
  expect_equivalent(dim(get_stats(scm.h5,per_chr = FALSE)),c(samples,4))
  
  expect_equivalent(dim(get_stats(scm.mem)),c(chr*samples,5))
  expect_equivalent(dim(get_stats(scm.mem,per_chr = FALSE)),c(samples,4))
  
})

test_that("get_region_summary", {
  
  expect_error(get_region_summary("not scMethrix"))
  expect_error(get_region_summary(scm.mem,group="not a group"))
  expect_error(get_region_summary(scm.mem,type="not a type"))
  expect_error(get_region_summary(scm.mem,how="not a how"))
  
  region <- GRanges(seqnames = c("chr1"), ranges = IRanges(1,10)) 
  expect_error(get_region_summary(scm.mem,region=region))
  
  region <- GRanges(seqnames = c("chr1","chr2","chr3","chr4"), ranges = IRanges(1,100000000)) 
  expect_equivalent(dim(get_region_summary(scm.mem,region=region)),c(4,9))
  #expect_warning(get_region_summary(scm.mem,n_chunks=1000,region=region))
  
})

test_that("mask_scMethrix", {
  
  expect_error(mask_scMethrix("not scMethrix"))
  expect_error(mask_scMethrix(scm.mem))
  expect_error(mask_scMethrix(scm.mem,n_threads=2))
  expect_error(mask_scMethrix(scm.mem,high_quantile=1,type="count"))
  expect_error(mask_scMethrix(scm.mem,high_quantile=5,type="coverage"))

  m <- mask_scMethrix(scm.mem,low_count=2)
  
  expect_equivalent(dim(m),c(100,4))
  expect_equivalent(dim(remove_uncovered(m)),c(75,4))
  
})
