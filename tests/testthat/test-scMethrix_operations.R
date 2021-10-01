test_that("get_metadata_stats", {
  
  expect_error(get_metadata_stats("not scMethrix"),msg.check.scm)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) { 
    expect_error(get_metadata_stats(scm="not scMethrix"))
    s <- get_metadata_stats(scm)
    expect_equal(dim(mcols(s)),c(n_cpg,5))
    
    s <- remove_assay(scm,assay="counts")
    s <- get_metadata_stats(s)
    expect_equal(dim(mcols(s)),c(n_cpg,4))
    
    stats <- mcols(s)
    rng <- 1:10
    
    expect_equal(rowMeans(score(s)[rng,],na.rm=TRUE),stats$mean_meth[rng])
    expect_equal(DelayedMatrixStats::rowMedians(score(s)[rng,],na.rm=TRUE),stats$median_meth[rng])
    expect_equal(DelayedMatrixStats::rowSds(score(s)[rng,],na.rm=TRUE),stats$sd_meth[rng])
    expect_equal(ncol(s)-rowCounts(score(s)[rng,],val=NA),stats$cells[rng])
  }))
})

test_that("remove_assay", {
  
  expect_error(remove_assay("not scMethrix"),msg.check.scm)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(remove_assay(scm, assay="not an assay"),msg.assay.match)
    expect_error(remove_assay(scm, assay="score"),"Score assay cannot be removed")
    
    plus1 <- transform_assay(scm,trans=function(x) x+1,assay="score",new_assay="plus1")
    expect_false(isTRUE(all.equal(assays(scm), assays(plus1))))
    plus1 <- remove_assay(plus1, assay="plus1")
    expect_true(all.equal(assays(scm), assays(plus1)))
  }))
})

test_that("merge_scMethrix", {
  
  expect_error(merge_scMethrix("not scMethrix"),msg.check.scm)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
   
    ### Assay tests
    ## Different assays
    expect_warning(merge_scMethrix(remove_assay(scm,assay="counts")[,1],scm[,2:ncol(scm)],by="col"),"Assay list not identical") 
    
    # Same assays
    if (!is_h5(scm)) { 
      expect_equal(merge_scMethrix(scm[1],scm[2:nrow(scm)],by="row"),scm)
      expect_equal(merge_scMethrix(scm[,1],scm[,2:ncol(scm)],by="col"),scm)
      
      #order doesn't matter
      expect_equal(merge_scMethrix(scm[1], scm[2:nrow(scm)],by="row"),
                   merge_scMethrix(scm[2:nrow(scm)], scm[1],by="row"))
      
      expect_equal(merge_scMethrix(scm[,1], scm[,2:ncol(scm)],by="col"),  #colbind is currently position dependant
                   merge_scMethrix(scm[,2:ncol(scm)], scm[,1],by="col"))
      
    } else {
      # Test is less stringent on HDF5-stored objects due since the filepath of assays cannot be equal
      expect_equal(as.matrix(score(merge_scMethrix(scm[1],scm[2:nrow(scm)],by="row"))),as.matrix(score(scm)))
      expect_equal(as.matrix(score(merge_scMethrix(scm[,1],scm[,2:ncol(scm)],by="col"))),as.matrix(score(scm)))
      
      #order doesn't matter
      expect_equal(as.matrix(score(merge_scMethrix(scm[1], scm[2:nrow(scm)],by="row"))),
                             as.matrix(score(merge_scMethrix(scm[2:nrow(scm)], scm[1],by="row"))))
      
      expect_equal(as.matrix(score(merge_scMethrix(scm[,1], scm[,2:ncol(scm)],by="col"))),
                   as.matrix(score(merge_scMethrix(scm[,2:ncol(scm)], scm[,1],by="col"))))
      
    }
    
    ### Metadata testing
    ## Column tests
    # Same samples
    expect_error(merge_scMethrix(scm,scm,by="col"), "You have the same samples in your datasets. ") 
    
    # Different regions
    expect_error(
      expect_warning(
        merge_scMethrix(scm[1,1],scm[2,2],by="col"),
        "Same metadata columns are present"
      ),"There are non-overlapping regions in your datasets.") 
    # Merging colData 
    scm1 <- scm[,1:2]; scm2 <- scm[,3:4]
    colData(scm1)$Sample <- colData(scm1)$ID <- 1; colData(scm2)$Sample <- 2
    scm12 = merge_scMethrix(scm1,scm2,by="col")
    expect_equal(colData(scm12)$Sample,c(rep(1,ncol(scm1)),rep(2,ncol(scm2))))
    expect_equal(colData(scm12)$ID,c(rep(1,ncol(scm1)),rep(NA,ncol(scm2))))
    # Merging mcols 
    mcols(scm1)$CpG <- mcols(scm1)$ID <- 1; mcols(scm2)$CpG <- 2
    expect_warning(scm12 <- merge_scMethrix(scm1,scm2,by="col"),"Same metadata columns are present")
    expect_true(setequal(names(mcols(scm12)),c("CpG.1","ID","CpG.2")))
    
    ## Row tests
    # Same regions
    expect_error(merge_scMethrix(scm,scm,by="row"),"There are overlapping regions in your datasets.") 
    # Different samples
    expect_error(
      expect_warning(
        merge_scMethrix(scm[1,1],scm[2,2],by="row"),
        "Same metadata columns are present"
      ),"You have different samples in your dataset.") 
    # Merging mcols
    scm1 <- scm[1:5]; scm2 <- scm[6:10]
    mcols(scm1)$CpG <- mcols(scm1)$ID <- 1; mcols(scm2)$CpG <- 2
    scm12 = merge_scMethrix(scm1,scm2,by="row")
    expect_equal(mcols(scm12)$CpG,c(rep(1,nrow(scm1)),rep(2,nrow(scm2))))
    expect_equal(mcols(scm12)$ID,c(rep(1,nrow(scm1)),rep(NA,nrow(scm2))))
    # Merging colData
    colData(scm1)$Sample <- colData(scm1)$ID <- 1; colData(scm2)$Sample <- 2
    expect_warning(scm12 <- merge_scMethrix(scm1,scm2,by="row"),"Same metadata columns are present")
    expect_true(setequal(names(colData(scm12)),c("Sample.1","ID","Sample.2")))
    
    # Check overall equality of merged object
    if (is_h5(scm)) {
      # Cannot test equality of objects directly, as the HDF5Matrix filepath will always differ between the two objects
      scm.mg1 <- scm
      scm.mg1 <- merge_scMethrix(scm.mg1[1],scm.mg1[2:nrow(scm)],by="row")
      scm.mg1 <- merge_scMethrix(scm.mg1[,1],scm.mg1[,2:ncol(scm)],by="col")
      
      expect_equal(metadata(scm.mg1),metadata(scm))
      expect_equal(colData(scm.mg1),colData(scm))
      expect_equal(mcols(scm.mg1),mcols(scm))
      expect_equal(rowRanges(scm.mg1),rowRanges(scm))
      expect_equal(as.matrix(score(scm.mg1)),as.matrix(score(scm)))
      
      # Make sure order doesn't matter
      scm.mg2 <- scm
      scm.mg2 <- merge_scMethrix(scm.mg2[2:nrow(scm)],scm.mg2[1],by="row")
      scm.mg2 <- merge_scMethrix(scm.mg2[,2:ncol(scm)],scm.mg2[,1],by="col")
      
      expect_equal(metadata(scm.mg1),metadata(scm.mg2))
      expect_equal(colData(scm.mg1),colData(scm.mg2))
      expect_equal(mcols(scm.mg1),mcols(scm.mg2))
      expect_equal(rowRanges(scm.mg1),rowRanges(scm.mg2))
      expect_equal(as.matrix(score(scm.mg1)),as.matrix(score(scm.mg2)))
      
    } else { 
      
      expect_equal(merge_scMethrix(scm[1],scm[2:nrow(scm)],by="row"),scm)
      expect_equal(merge_scMethrix(scm[,1],scm[,2:ncol(scm)],by="col"),scm)
      
      # Make sure order doesn't matter
      expect_equal(merge_scMethrix(scm[1], scm[2:nrow(scm)],by="row"),
                   merge_scMethrix(scm[2:nrow(scm)], scm[1],by="row"))
      expect_equal(merge_scMethrix(scm[,1],scm[,2:ncol(scm)],by="col"),
                   merge_scMethrix(scm[,2:ncol(scm)],scm[,1],by="col"))
    }
    
  }))
})

test_that("convert_HDF5_scMethrix", {

  expect_error(convert_HDF5_scMethrix("not scMethrix"),msg.check.scm)
  
  expect_true(is_h5(scm.h5))
  expect_equal(class(get_matrix(scm.h5))[1],"HDF5Matrix")
  
  scm <- convert_HDF5_scMethrix(scm.h5)
  
  expect_false(is_h5(scm))
  expect_equal(as.numeric(score(scm)),as.numeric(score(scm.h5)))
  expect_equal(class(get_matrix(scm))[1],"matrix") 
  rm(scm)
  
})

test_that("convert_scMethrix", {
  
  expect_error(convert_scMethrix("not scMethrix"),msg.check.scm)
  
  expect_false(is_h5(scm.mem))
  expect_equal(class(get_matrix(scm.mem))[1],"matrix") 
  
  scm <- convert_scMethrix(scm.mem)
  
  expect_true(is_h5(scm))
  expect_equal(as.numeric(score(scm)),as.numeric(score(scm.mem)))
  expect_equal(class(get_matrix(scm))[1],"HDF5Matrix")
  rm(scm)
  
})

test_that("subset_scMethrix", {
  
  expect_error(subset_scMethrix("not scMethrix"),msg.check.scm)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_warning(subset_scMethrix(scm))
    
    samples <- c("C1","C3")
    contigs <- c("chr1")
    regions <- GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(1,100000000)) 
    
    # Subset by include
    s <- subset_scMethrix(scm, samples = samples, by="include")
    expect_equal(dim(s),c(n_cpg,length(samples)))
    expect_equal(samples,colData(s)@rownames)
   
    s <- subset_scMethrix(scm, contigs = contigs, by="include")
    expect_equal(dim(s),c(147,n_samples))
    expect_equal(contigs,as.character(seqnames(s)@values))
    
    s <- subset_scMethrix(scm, regions = regions, by="include")
    expect_equal(dim(s),c(134,n_samples))
    expect_equal(length(findOverlaps(regions,rowRanges(s))),length(rowRanges(s)))
      
    s <- subset_scMethrix(scm, samples = samples, contigs = contigs, regions = regions, by="include")
    expect_equal(dim(s),c(67,length(samples)))
    
    # Subset by exclude
    s <- subset_scMethrix(scm, samples = samples, by = "exclude")
    expect_equal(dim(s),c(n_cpg,length(samples)))
    expect_equal(length(intersect(rownames(colData(s)),samples)),0)
    
    s <- subset_scMethrix(scm, contigs = contigs, by = "exclude")
    expect_equal(dim(s),c(139,n_samples))
    expect_equal(length(intersect(contigs,as.character(seqnames(s)@values))),0)
    
    s <- subset_scMethrix(scm, regions = regions, by = "exclude")
    expect_equal(dim(s),c(152,n_samples))
    expect_equal(length(findOverlaps(regions,rowRanges(s))),0)
    
    s <- subset_scMethrix(scm, samples = samples, contigs = contigs, regions = regions, by = "exclude")
    expect_equal(dim(s),c(72,length(samples)))
    
  }))
})

test_that("get_matrix", {
 
  expect_error(get_matrix("not scMethrix"),msg.check.scm)
  
    s <- get_matrix(scm.h5)
    expect_equal(dim(s),c(n_cpg,4))  
    expect_equal(class(s)[1],"HDF5Matrix")
    
    s <- get_matrix(scm.mem)
    expect_equal(dim(s),c(n_cpg,4))  
    expect_equal(class(s)[1],"matrix")
    
    invisible(lapply(list(scm.mem,scm.h5), function(scm) {
      expect_warning(get_matrix(scm,add_loci=FALSE, in_granges = TRUE))
      
      s <- get_matrix(scm=scm,add_loci=TRUE)
      expect_equal(dim(s),c(n_cpg,7))  
      expect_equal(class(s)[1],"data.table")
      
      s <- get_matrix(scm,add_loci=TRUE, in_granges = TRUE)
     # expect_equal(seqnames(m)@lengths,c(10,8))
      expect_equal(dim(mcols(s)),c(n_cpg,4))
      expect_equal(class(s)[1],"GRanges")
      
      s <- get_matrix(scm,order_by_sd = TRUE)
      expect_false(is.unsorted(rev(rowSds(s,na.rm=TRUE)),na.rm=TRUE))
  }))
})


test_that("remove_uncovered", {
  
  expect_error(remove_uncovered("not scMethrix"),msg.check.scm)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    samples = c("C1","C2")
    s <- subset_scMethrix(scm,samples=samples,by = "include")
    expect_equal(dim(s),c(n_cpg,length(samples)))
    non_na_rows = nrow(score(s)[rowSums(is.na(score(s))) != ncol(s), ])
    expect_equal(dim(remove_uncovered(s)),c(non_na_rows,length(samples)))
  }))
})

test_that("get_stats", {
  
  expect_error(get_stats("not scMethrix"),msg.check.scm)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    chr <- length(seqlengths(rowRanges(scm)))
    samples <- nrow(colData(scm))
    expect_equal(dim(get_stats(scm)),c(chr*samples,5))
    expect_equal(dim(get_stats(scm,per_chr = FALSE)),c(samples,4))
    
    stats <- get_stats(scm,per_chr=FALSE)
    expect_equal(mean(score(scm)[,"C1"],na.rm=TRUE), as.double(stats[Sample_Name == "C1","mean_meth"]))
    expect_equal(median(score(scm)[,"C1"],na.rm=TRUE), as.double(stats[Sample_Name == "C1","median_meth"]))
    expect_equal(sd(score(scm)[,"C1"],na.rm=TRUE), as.double(stats[Sample_Name == "C1","sd_meth"]))
    
    scm <- subset_scMethrix(scm,contigs="chr1")
    stats <- get_stats(scm,per_chr=TRUE)
    expect_equal(mean(score(scm)[,"C1"],na.rm=TRUE), as.double(stats[Sample_Name == "C1","mean_meth"]))
    expect_equal(median(score(scm)[,"C1"],na.rm=TRUE), as.double(stats[Sample_Name == "C1","median_meth"]))
    expect_equal(sd(score(scm)[,"C1"],na.rm=TRUE), as.double(stats[Sample_Name == "C1","sd_meth"]))
    
  }))
})

test_that("get_region_summary", {
  
  expect_error(get_region_summary("not scMethrix"),msg.check.scm)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(get_region_summary(scm,group="not a group"))
    expect_error(get_region_summary(scm,type="not a type"))
    expect_error(get_region_summary(scm,how="not a how"))
    
    region <- GRanges(seqnames = c("chr1"), ranges = IRanges(1,10)) 
    expect_error(get_region_summary(scm,region=region))
    
    region <- GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(1,100000000)) 
    expect_equal(dim(get_region_summary(scm,region=region)),c(2,9))
  #expect_warning(get_region_summary(scm.mem,n_chunks=1000,region=region))
  }))
})

test_that("mask_by_coverage", {
  expect_error(mask_by_coverage("not scMethrix"),msg.check.scm)
  expect_error(mask_by_coverage(scm.mem,assay="not an assay"),"Invalid assay")
  expect_error(mask_by_coverage(remove_assay(scm.mem,assay="counts")))
  expect_error(mask_by_coverage(scm.mem,n_threads=2))
  expect_error(mask_by_coverage(scm.mem,low_threshold=-1),"low_threshold")
  expect_error(mask_by_coverage(scm.mem,low_threshold="not numeric"),"low_threshold")
  expect_error(mask_by_coverage(scm.mem,avg_threshold=-1,"avg_threshold"))
  expect_error(mask_by_coverage(scm.mem,avg_threshold="not numeric"),"avg_threshold")
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {

    scm <- mask_by_coverage(scm,low_threshold=2,avg_threshold=NULL)
    expect_equal(dim(scm),c(n_cpg,n_samples))
    expect_equal(dim(remove_uncovered(scm)),c(length(which(rowSums(counts(scm),na.rm=TRUE) >= 2)),n_samples))

    scm <- mask_by_coverage(scm,low_threshold=NULL,avg_threshold=1) #Removes all rows with a sample with a coverage of 2
    expect_equal(dim(scm),c(n_cpg,n_samples))
    expect_equal(dim(remove_uncovered(scm)),c(length(which(rowMeans(counts(scm),na.rm=TRUE) == 1)),n_samples))

  }))
})

test_that("mask_by_sample", {
  expect_error(mask_by_sample("not scMethrix"),msg.check.scm)
  expect_error(mask_by_sample(scm.mem,assay="not an assay"),"Invalid assay")
  expect_error(mask_by_sample(scm.mem,n_threads=2))
  expect_error(mask_by_sample(scm.mem,low_threshold=2,prop_threshold=1))
  expect_error(mask_by_sample(scm.mem,low_threshold=-1),"low_threshold")
  expect_error(mask_by_sample(scm.mem,low_threshold="not numeric"),"low_threshold")
  expect_error(mask_by_sample(scm.mem,prop_threshold=-1,"prop_threshold"))
  expect_error(mask_by_sample(scm.mem,prop_threshold=2,"prop_threshold"))
  expect_error(mask_by_sample(scm.mem,prop_threshold="not numeric"),"prop_threshold")
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    scm <- mask_by_sample(scm,low_threshold=2)
    expect_equal(dim(scm),c(n_cpg,n_samples))
    expect_equal(dim(remove_uncovered(scm)),c(length(which(ncol(scm) - rowCounts(score(scm),value=NA) >= 2)),n_samples))
    
    scm <- mask_by_sample(scm,low_threshold=NULL,prop_threshold=0.25) # Since 0.25 of 4 sample is 1, same result as low_threshold = 1
    expect_equal(dim(scm),c(n_cpg,n_samples))
    expect_equal(dim(remove_uncovered(scm)),c(length(which(ncol(scm) - rowCounts(score(scm),value=NA) > 1)),n_samples))
    
  }))
})

test_that("mask_by_variance", {
  expect_error(mask_by_variance("not scMethrix"),msg.check.scm)
  expect_error(mask_by_variance(scm.mem,assay="not an assay"),"Invalid assay")
  expect_error(mask_by_variance(scm.mem,n_threads=2))
  expect_error(mask_by_variance(scm.mem,low_threshold=2,"low_threshold must be between 0 and 1"))
  expect_error(mask_by_variance(scm.mem,low_threshold=-1,"low_threshold must be between 0 and 1"))
  expect_error(mask_by_variance(scm.mem,low_threshold="not numeric"),"low_threshold must be between 0 and 1")
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    #Since only 4 samples, only rows with 1 unique score value are masked
    uniq <- apply(score(scm),1,function(x) {
      uniq <- unique(x)
      return(length(uniq[!is.na(uniq)]))
    })
    
    scm <- mask_by_variance(scm,low_threshold=0.05)

    expect_equal(dim(scm),c(n_cpg,n_samples))
    expect_equal(dim(remove_uncovered(scm)),c(sum(uniq != 1),n_samples))
  }))
})