test_that("convert_HDF5_methrix", {

  expect_error(convert_HDF5_methrix("not scMethrix"))
  
  expect_true(is_h5(scm.h5))
  expect_equivalent(class(get_matrix(scm.h5))[1],"HDF5Matrix")
  
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
  contigs <- c("chr1")
  regions <- GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(1,8)) 
  
  # Subset H5 by include
  
  s <- subset_scMethrix(scm.h5, samples = samples, by="include")
  expect_equivalent(dim(s),c(18,2))
  expect_equivalent(samples,colData(s)$colData)

  s <- subset_scMethrix(scm.h5, contigs = contigs, by="include")
  expect_equivalent(dim(s),c(10,4))
  expect_equivalent(contigs,as.character(seqnames(s)@values))
  
  s <- subset_scMethrix(scm.h5, regions = regions, by="include")
  expect_equivalent(dim(s),c(6,4))

  s <- subset_scMethrix(scm.h5, samples = samples, contigs = contigs, regions = regions, by="include")
  expect_equivalent(dim(s),c(4,2))

  # Subset mem by include
  
  s <- subset_scMethrix(scm.mem, samples = samples, by="include")
  expect_equivalent(dim(s),c(18,2))

  s <- subset_scMethrix(scm.mem, contigs = contigs, by="include")
  expect_equivalent(dim(s),c(10,4))
  
  s <- subset_scMethrix(scm.mem, regions = regions, by="include")
  expect_equivalent(dim(s),c(6,4))
  
  s <- subset_scMethrix(scm.mem, samples = samples, contigs = contigs, regions = regions, by="include")
  expect_equivalent(dim(s),c(4,2))
  
  # Subset H5 by exclude
  
  s <- subset_scMethrix(scm.h5, samples = samples, by = "exclude")
  expect_equivalent(dim(s),c(18,2))
  expect_equivalent(length(intersect(colData(s)$colData,samples)),0)
  
  s <- subset_scMethrix(scm.h5, contigs = contigs, by = "exclude")
  expect_equivalent(dim(s),c(8,4))
  #expect_equivalent(contigs,as.character(seqnames(s)@values))
  
  s <- subset_scMethrix(scm.h5, regions = regions, by = "exclude")
  expect_equivalent(dim(s),c(12,4))
  
  s <- subset_scMethrix(scm.h5, samples = samples, contigs = contigs, regions = regions, by = "exclude")
  expect_equivalent(dim(s),c(6,2))
  
  # Subset mem by exclude
  
  s <- subset_scMethrix(scm.mem, samples = samples, by = "exclude")
  expect_equivalent(dim(s),c(18,2))
  
  s <- subset_scMethrix(scm.mem, contigs = contigs, by = "exclude")
  expect_equivalent(dim(s),c(8,4))
  
  s <- subset_scMethrix(scm.mem, regions = regions, by = "exclude")
  expect_equivalent(dim(s),c(12,4))
  
  s <- subset_scMethrix(scm.mem, samples = samples, contigs = contigs, regions = regions, by = "exclude")
  expect_equivalent(dim(s),c(6,2))
  
  
})

test_that("get_matrix", {

  expect_error(get_matrix("not scMethrix"))
  expect_warning(get_matrix(scm.h5,add_loci=FALSE, in_granges = TRUE))

  m <- get_matrix(scm.h5)
  expect_equivalent(dim(m),c(18,4))  
  expect_equivalent(class(m)[1],"HDF5Matrix")
  
  m <- get_matrix(m=scm.h5,add_loci=TRUE)
  expect_equivalent(dim(m),c(18,7))  
  expect_equivalent(class(m)[1],"data.table")
  
  m <- get_matrix(scm.h5,add_loci=TRUE, in_granges = TRUE)
  expect_equivalent(seqnames(m)@lengths,c(10,8))
  expect_equivalent(dim(mcols(m)),c(18,4))
  expect_equivalent(class(m)[1],"GRanges")

  m <- get_matrix(scm.mem)
  expect_equivalent(dim(m),c(18,4))  
  expect_equivalent(class(m)[1],"matrix")
  
  m <- get_matrix(scm.mem,add_loci=TRUE)
  expect_equivalent(dim(m),c(18,7))  
  expect_equivalent(class(m)[1],"data.table")
  
  m <- get_matrix(scm.mem,add_loci=TRUE, in_granges = TRUE)
  expect_equivalent(seqnames(m)@lengths,c(10,8))
  expect_equivalent(dim(mcols(m)),c(18,4))
  expect_equivalent(class(m)[1],"GRanges")
  
})


test_that("remove_uncovered", {
  
  expect_error(get_matrix("not scMethrix"))

  h5 <- subset_scMethrix(scm.h5,samples="df1",by = "include")
  expect_equivalent(dim(h5),c(18,1))
  expect_equivalent(dim(remove_uncovered(h5)),c(5,1))

  mem <- subset_scMethrix(scm.mem,samples="df1",by = "include")
  expect_equivalent(dim(h5),c(18,1))
  expect_equivalent(dim(remove_uncovered(mem)),c(5,1))
  
  expect_equivalent(h5,mem)
  
  rm(h5,mem)
  
})

test_that("get_stats", {
  
  expect_error(get_matrix("not scMethrix"))

  expect_equivalent(dim(get_stats(scm.h5)),c(8,5))
  expect_equivalent(dim(get_stats(scm.h5,per_chr = FALSE)),c(4,4))
  
  expect_equivalent(dim(get_stats(scm.mem)),c(8,5))
  expect_equivalent(dim(get_stats(scm.mem,per_chr = FALSE)),c(4,4))
  
})

test_that("get_region_summary", {
  
  expect_error(get_region_summary("not scMethrix"))
  expect_error(get_region_summary(scm.mem,group="not a group"))
  expect_warning(get_region_summary(scm.mem,n_chunks=1000))
  expect_error(get_region_summary(scm.mem,type="not a type"))
  expect_error(get_region_summary(scm.mem,how="not a how"))
  
  region <- GRanges(seqnames = c("chr3"), ranges = IRanges(1,10)) 
  expect_error(get_region_summary(scm.mem,region=region))
  
  region <- GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(1,8)) 
  expect_equivalent(dim(get_region_summary(scm.mem,region=region)),c(2,9))

})

test_that("mask_methrix", {
  
  expect_error(mask_methrix("not scMethrix"))
  expect_error(mask_methrix(scm.mem))
  expect_error(mask_methrix(scm.mem,n_threads=2))
  expect_error(mask_methrix(scm.mem,high_quantile=1,type="count"))
  expect_error(mask_methrix(scm.mem,high_quantile=5,type="coverage"))
  
  expect_error(get_region_summary(scm.mem,type="not a type"))
  expect_error(get_region_summary(scm.mem,how="not a how"))
  
  expect_equivalent(dim(mask_methrix(scm.mem,low_count=10)),c(18,4))
  
})
