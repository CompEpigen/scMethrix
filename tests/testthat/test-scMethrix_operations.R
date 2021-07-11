test_that("get_metadata_stats", {
  invisible(lapply(list(scm.mem,scm.h5), function(scm) { 
    expect_error(get_metadata_stats(scm="not scMethrix"))
    s <- get_metadata_stats(scm)
    expect_equivalent(dim(mcols(s)),c(n_cpg,5))
    
    s <- remove_assay(scm,assay="counts")
    s <- get_metadata_stats(s)
    expect_equivalent(dim(mcols(s)),c(n_cpg,4))
  }))
})

test_that("remove_assay", {
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(remove_assay(scm="not scMethrix"))
    expect_error(remove_assay(scm, assay="not an assay"))
    expect_error(remove_assay(scm, assay="score"))
    plus1 <- transform_assay(scm,trans=function(x) x+1,assay="score",new_assay="plus1")
    expect_false(isTRUE(all.equal(assays(scm), assays(plus1))))
    plus1 <- remove_assay(plus1, assay="plus1")
    expect_true(all.equal(assays(scm), assays(plus1)))
  }))
})

test_that("merge_scMethrix", {
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(merge_scMethrix(scm1=scm,scm2="not scMethrix"))
    
    expect_error(merge_scMethrix(scm,scm,by="col")) #same samples
    expect_error(merge_scMethrix(scm[1,1],scm[2,2],by="col")) #different regions
    
    expect_error(merge_scMethrix(scm,scm,by="row")) #same regions
    expect_error(merge_scMethrix(scm[1,1],scm[2,2],by="row")) #different samples
    
    s <- scm
    assays(s) <- assays(s)["score"]
    expect_warning(merge_scMethrix(s[,1],scm[,2],by="col")) #different assays
    
    expect_equivalent(merge_scMethrix(scm[1],scm[2:nrow(scm)],by="row"),scm)
    expect_equivalent(merge_scMethrix(scm[,1],scm[,2:ncol(scm)],by="col"),scm)
  }))
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
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_warning(subset_scMethrix(scm))
    
    samples <- c("C1","C3")
    contigs <- c("chr1")
    regions <- GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(1,100000000)) 
    
    # Subset by include
    s <- subset_scMethrix(scm, samples = samples, by="include")
    expect_equivalent(dim(s),c(n_cpg,length(samples)))
    expect_equivalent(samples,colData(s)@rownames)
   
    s <- subset_scMethrix(scm, contigs = contigs, by="include")
    expect_equivalent(dim(s),c(147,n_samples))
    expect_equivalent(contigs,as.character(seqnames(s)@values))
    
    s <- subset_scMethrix(scm, regions = regions, by="include")
    expect_equivalent(dim(s),c(134,n_samples))
  
    s <- subset_scMethrix(scm, samples = samples, contigs = contigs, regions = regions, by="include")
    expect_equivalent(dim(s),c(67,length(samples)))
    
    # Subset by exclude
    s <- subset_scMethrix(scm, samples = samples, by = "exclude")
    expect_equivalent(dim(s),c(n_cpg,length(samples)))
    expect_equivalent(length(intersect(colData(s)$colData,samples)),0)
    
    s <- subset_scMethrix(scm, contigs = contigs, by = "exclude")
    expect_equivalent(dim(s),c(139,n_samples))
    #expect_equivalent(contigs,as.character(seqnames(s)@values))
    
    s <- subset_scMethrix(scm, regions = regions, by = "exclude")
    expect_equivalent(dim(s),c(152,n_samples))
    
    s <- subset_scMethrix(scm, samples = samples, contigs = contigs, regions = regions, by = "exclude")
    expect_equivalent(dim(s),c(72,length(samples)))
  }))
})

test_that("get_matrix", {
 
    expect_error(get_matrix("not scMethrix"))
  
    s <- get_matrix(scm.h5)
    expect_equivalent(dim(s),c(n_cpg,4))  
    expect_equivalent(class(s)[1],"HDF5Matrix")
    
    s <- get_matrix(scm.mem)
    expect_equivalent(dim(s),c(n_cpg,4))  
    expect_equivalent(class(s)[1],"matrix")
    
    invisible(lapply(list(scm.mem,scm.h5), function(scm) {
      expect_warning(get_matrix(scm,add_loci=FALSE, in_granges = TRUE))
      
      s <- get_matrix(scm=scm,add_loci=TRUE)
      expect_equivalent(dim(s),c(n_cpg,7))  
      expect_equivalent(class(s)[1],"data.table")
      
      s <- get_matrix(scm,add_loci=TRUE, in_granges = TRUE)
     # expect_equivalent(seqnames(m)@lengths,c(10,8))
      expect_equivalent(dim(mcols(s)),c(n_cpg,4))
      expect_equivalent(class(s)[1],"GRanges")
      
      s <- get_matrix(scm,order_by_sd = TRUE)
      expect_false(is.unsorted(rev(rowSds(s,na.rm=TRUE)),na.rm=TRUE))
  }))
})


test_that("remove_uncovered", {
  
  expect_error(get_matrix("not scMethrix"))
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    samples = c("C1","C2")
    s <- subset_scMethrix(scm,samples=samples,by = "include")
    expect_equivalent(dim(s),c(n_cpg,length(samples)))
    expect_equivalent(dim(remove_uncovered(s)),c(235,2))
  }))
})

test_that("get_stats", {
  
  expect_error(get_matrix("not scMethrix"))
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    chr <- length(seqlengths(rowRanges(scm)))
    samples <- nrow(colData(scm))
    expect_equivalent(dim(get_stats(scm)),c(chr*samples,5))
    expect_equivalent(dim(get_stats(scm,per_chr = FALSE)),c(samples,4))
  }))
})

test_that("get_region_summary", {
  
  expect_error(get_region_summary("not scMethrix"))
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(get_region_summary(scm,group="not a group"))
    expect_error(get_region_summary(scm,type="not a type"))
    expect_error(get_region_summary(scm,how="not a how"))
    
    region <- GRanges(seqnames = c("chr1"), ranges = IRanges(1,10)) 
    expect_error(get_region_summary(scm,region=region))
    
    region <- GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(1,100000000)) 
    expect_equivalent(dim(get_region_summary(scm,region=region)),c(2,9))
  #expect_warning(get_region_summary(scm.mem,n_chunks=1000,region=region))
  }))
})

test_that("mask_scMethrix", {
  
  expect_error(mask_scMethrix("not scMethrix"))
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    #expect_error(mask_scMethrix(scm,n_threads=2))
    #expect_error(mask_scMethrix(scm,max_avg_count=1,type="cells"))
    
    m <- mask_scMethrix(scm,low_total_count=2, type="counts")
    expect_equivalent(dim(m),c(n_cpg,n_samples))
    expect_equivalent(dim(score(m)),dim(counts(m)))
    expect_equivalent(dim(remove_uncovered(m)),c(232,n_samples))
    
    m <- mask_scMethrix(scm,max_avg_count=1,type="counts")
    expect_equivalent(dim(m),c(n_cpg,n_samples))
    expect_equivalent(dim(score(m)),dim(counts(m)))
    expect_equivalent(dim(remove_uncovered(m)),c(170,n_samples))
    
    m <- mask_scMethrix(scm,low_total_count=2,type="cells")
    expect_equivalent(dim(m),c(n_cpg,n_samples))
    expect_equivalent(dim(remove_uncovered(m)),c(219,n_samples))
    
  }))
})
