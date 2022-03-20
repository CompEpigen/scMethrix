test_that("save_scMethrix", {
  
  expect_error(save_scMethrix("not scMethrix"),"A valid SummarizedExperiment-derived object needs to be supplied.")
  expect_warning(save_scMethrix(scm.h5),"No dest specified")

  temp_dir <- tempfile("save_")
  scm.saved <- save_scMethrix(scm.h5,dest = temp_dir)
  expect_s4_class(scm.saved,"scMethrix")
  
  # Simulate the user pressing no on the replace folder dialog
  local({
    local_mock(menu = function(choices,title=NULL) 2)
    
    test_that("save_scMethrix | Replace: No", {
      expect_message(save_scMethrix(scm.h5,dest = temp_dir), "Saving aborted.")
    })
  })
  
  # Simulate the user pressing yes on the replace folder dialog
  local({
    local_mock(menu = function(choices,title=NULL) 1)
    
    test_that("save_scMethrix | Replace: Yes", {
      bad.file = paste0(temp_dir,"/bad_file.txt")
      file.create(bad.file)
      expect_true(file.exists(bad.file))
      scm.repl <- save_scMethrix(scm.h5,dest = temp_dir)
      expect_s4_class(scm.repl,"scMethrix")
      expect_false(file.exists(bad.file))
    })
  })
  
  bad.file = paste0(temp_dir,"/bad_file.txt")
  file.create(bad.file)
  expect_true(file.exists(bad.file))
  scm.saved <- save_scMethrix(scm.h5,dest = temp_dir,replace=TRUE)
  expect_s4_class(scm.saved,"scMethrix")
  expect_false(file.exists(bad.file))
  
  scm.quick <- save_scMethrix(scm.saved,quick=TRUE) #TODO: Should try and capture output to see what's happening
  expect_s4_class(scm.quick,"scMethrix")
  
  expect_warning(save_scMethrix(scm.saved,dest = temp_dir,quick=TRUE),"dest is not used ")
  
})

test_that("get_rowdata_stats", {
  
  fmt <- function(x) round(as.numeric(x),4)
  
  expect_error(get_rowdata_stats("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) { 
    
    cols <- ncol(rowData(scm))
    s <- get_rowdata_stats(scm)
    expect_equal(dim(rowData(s)),c(n_cpg,cols+4))
    stats <- rowData(s)
    
    # Check stats
    expect_equal(fmt(DelayedMatrixStats::rowMeans2(score(s),na.rm=TRUE)),fmt(stats$Mean))
    expect_equal(fmt(ncol(scm)-DelayedMatrixStats::rowCounts(score(s),val=NA)),fmt(stats$Cells))
    expect_equal(fmt((ncol(scm)-DelayedMatrixStats::rowCounts(score(s),val=NA))/ncol(scm)),fmt(stats$Sparsity))
    
    exp_sd <- DelayedMatrixStats::rowSds(score(s),na.rm=TRUE)
    exp_sd[is.na(exp_sd)] <- 0 # Since single sample rows will give SD as zero
    expect_equal(fmt(exp_sd),fmt(stats$SD))
    
    # Check for stat subset
    rd <- rowData(get_rowdata_stats(scm,stats=c("Mean")))
    expect_equal(dim(rd),c(n_cpg,cols+1))
    expect_true(colnames(rd)[cols+1] == "Mean")
    rd <- rowData(get_rowdata_stats(scm,stats=c("SD")))
    expect_true(colnames(rd)[cols+1] == "SD")
    rd <- rowData(get_rowdata_stats(scm,stats=c("Cells")))
    expect_true(colnames(rd)[cols+1] == "Cells")
    rd <- rowData(get_rowdata_stats(scm,stats=c("Sparsity")))
    expect_true(colnames(rd)[cols+1] == "Sparsity")
    
    rd <- rowData(get_rowdata_stats(scm,stat=c("Mean","SD")))
    expect_equal(dim(rd),c(n_cpg,cols+2))
    expect_true(any(colnames(rd) == "Mean"))
    expect_true(any(colnames(rd) =="SD"))
    
  }))
})

test_that("get_coldata_stats", {
  
  fmt <- function(x) round(as.numeric(x),4)
  
  expect_error(get_coldata_stats("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) { 
    
    cols <- ncol(rowData(scm))
    s <- get_coldata_stats(scm)
    expect_equal(dim(colData(s)),c(n_samples,cols+4))
    stats <- colData(s)
    
    # Check stats
    expect_equal(fmt(DelayedMatrixStats::colMeans2(score(s),na.rm=TRUE)),fmt(stats$Mean))
    #expect_equal(DelayedMatrixStats::rowMedians(score(s)[rng,],na.rm=TRUE),stats$median[rng])
    expect_equal(fmt(DelayedMatrixStats::colSds(score(s),na.rm=TRUE)),fmt(stats$SD))
    expect_equal(fmt(nrow(scm)-DelayedMatrixStats::colCounts(score(s),val=NA)),fmt(stats$CpGs))
    expect_equal(fmt((nrow(scm)-DelayedMatrixStats::colCounts(score(s),val=NA))/nrow(scm)),fmt(stats$Sparsity))
    
    # Check for stat subset
    cd <- colData(get_coldata_stats(scm,stats=c("Mean")))
    expect_equal(dim(cd),c(n_samples,cols+1))
    expect_true(colnames(cd)[cols+1] == "Mean")
    cd <- colData(get_coldata_stats(scm,stats=c("SD")))
    expect_true(colnames(cd)[cols+1] =="SD")
    cd <- colData(get_coldata_stats(scm,stats=c("CpGs")))
    expect_true(colnames(cd)[cols+1] == "CpGs")
    cd <- colData(get_coldata_stats(scm,stats=c("Sparsity")))
    expect_true(colnames(cd)[cols+1] == "Sparsity")
    
    cd <- colData(get_coldata_stats(scm,stat=c("Mean","SD")))
    expect_equal(dim(cd),c(n_samples,cols+2))
    expect_true(any(colnames(cd) == "Mean"))
    expect_true(any(colnames(cd) == "SD"))
  }))
})

test_that("remove_assay", {
  
  expect_error(remove_assay("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(remove_assay(scm, assay="not an assay"),msg.validateAssay)
    expect_error(remove_assay(scm, assay="score"),"Score assay cannot be removed")
    
    plus1 <- transform_assay(scm,trans=function(x) x+1,assay="score",new_assay="plus1")
    expect_false(isTRUE(all.equal(assays(scm), assays(plus1))))
    plus1 <- remove_assay(plus1, assay="plus1")
    expect_true(all.equal(assays(scm), assays(plus1)))
  }))
})

test_that("merge_scMethrix", {
  
  expect_error(merge_scMethrix("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
   
    ### Assay tests
    ## Different assays
    expect_warning(merge_scMethrix(remove_assay(scm,assay="counts")[,1],scm[,2:ncol(scm)],by="col"),"Assay list not identical") 
    
    # Same assays
    identical.scm(merge_scMethrix(scm[1],scm[2:nrow(scm)],by="row"),scm)
    identical.scm(merge_scMethrix(scm[,1],scm[,2:ncol(scm)],by="col"),scm)
    
    # Order doesn't matter
    identical.scm(merge_scMethrix(scm[1], scm[2:nrow(scm)],by="row"),
                  merge_scMethrix(scm[2:nrow(scm)], scm[1],by="row"))
    identical.scm(merge_scMethrix(scm[,1], scm[,2:ncol(scm)],by="col"),  #cbindlist is currently position dependant
                  merge_scMethrix(scm[,2:ncol(scm)], scm[,1],by="col"))
    

    # Check H5 separately
    if (is_h5(scm)) { 
      expect_true(is_h5(merge_scMethrix(scm[1],scm[2:nrow(scm)],by="row")))
    } else {
      expect_false(is_h5(merge_scMethrix(scm[1],scm[2:nrow(scm)],by="row")))
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
    scm1 <- scm[,1:2]
    scm2 <- scm[,3:4]
    colData(scm1)$Group <- colData(scm1)$ID <- 1; colData(scm2)$Group <- 2
    scm12 = merge_scMethrix(scm1,scm2,by="col")
    expect_equal(colData(scm12)$Group,c(rep(1,ncol(scm1)),rep(2,ncol(scm2))))
    expect_equal(colData(scm12)$ID,c(rep(1,ncol(scm1)),rep(NA,ncol(scm2))))
    # Merging rowData 
    rowData(scm1)$CpG <- rowData(scm1)$ID <- 1; rowData(scm2)$CpG <- 2
    expect_warning(scm12 <- merge_scMethrix(scm1,scm2,by="col"),"Same metadata columns are present")
    expect_true(setequal(names(rowData(scm12)),c("CpG.1","ID","CpG.2")))
    
    ## Row tests
    # Same regions
    expect_error(merge_scMethrix(scm,scm,by="row"),"There are overlapping regions in your datasets.") 
    # Different samples
    # expect_error(
    #   expect_warning(
    #     merge_scMethrix(scm[1,1],scm[2,2],by="row"),
    #     "Same metadata columns are present"
    #   ),"You have different samples in your dataset.") 
    # Merging mcols
    scm1 <- scm[1:5]
    scm2 <- scm[6:10]
    rowData(scm1)$CpG <- rowData(scm1)$ID <- 1; rowData(scm2)$CpG <- 2
    scm12 = merge_scMethrix(scm1,scm2,by="row")
    expect_equal(rowData(scm12)$CpG,c(rep(1,nrow(scm1)),rep(2,nrow(scm2))))
    expect_equal(rowData(scm12)$ID,c(rep(1,nrow(scm1)),rep(NA,nrow(scm2))))
    # Merging colData 
    colData(scm1)$Group <- colData(scm1)$ID <- 1; colData(scm2)$Group <- 2
    expect_warning(scm12 <- merge_scMethrix(scm1,scm2,by="row"),"Same metadata columns are present")
    expect_true(setequal(names(colData(scm12)),c("Group.1","ID","Group.2")))
    
  }))
})

test_that("convert_scMethrix", {
  
  expect_error(convert_scMethrix("not scMethrix"),msg.validateExp)
  
  # convert memory to hdf5
  
  expect_false(is_h5(scm.mem))
  expect_equal(class(get_matrix(scm.mem))[1],"matrix") 
  
  scm <- convert_scMethrix(scm.mem, h5_dir = paste0(tempdir(),"/h5"))
  
  expect_true(is_h5(scm))
  
  for (name in assayNames(scm)) {
    expect_equal(as.numeric(get_matrix(scm,name)),as.numeric(get_matrix(scm.h5,name)))
    expect_true("HDF5Matrix" %in% class(get_matrix(scm,name))) 
  }
  
  # convert hdf5 to memory
  
  expect_true(is_h5(scm.h5))
  expect_equal(class(get_matrix(scm.h5))[1],"HDF5Matrix")
  
  scm <- convert_scMethrix(scm.h5)
  
  expect_false(is_h5(scm))
  
  for (name in assayNames(scm)) {
    expect_equal(as.numeric(get_matrix(scm,name)),as.numeric(get_matrix(scm.h5,name)))
    expect_true("matrix" %in% class(get_matrix(scm,name))) 
  }
  
  # convert HDF5, but do nothing because it's already in type
  
  expect_true(is_h5(scm.h5))
  expect_equal(class(get_matrix(scm.h5))[1],"HDF5Matrix")
  
  scm <- convert_scMethrix(scm.h5,type="HDF5")
  
  expect_true(is_h5(scm))
  
  for (name in assayNames(scm)) {
    expect_equal(as.numeric(get_matrix(scm,name)),as.numeric(get_matrix(scm.h5,name)))
    expect_true("HDF5Matrix" %in% class(get_matrix(scm,name))) 
  }
  
  # convert memory, but do nothing because it's already in type
  
  expect_false(is_h5(scm.mem))
  expect_equal(class(get_matrix(scm.mem))[1],"matrix") 
  
  scm <- convert_scMethrix(scm.mem, type = "memory", h5_dir = paste0(tempdir(),"/h5"))
  
  expect_false(is_h5(scm))
  
  for (name in assayNames(scm)) {
    expect_equal(as.numeric(get_matrix(scm,name)),as.numeric(get_matrix(scm.h5,name)))
    expect_true("matrix" %in% class(get_matrix(scm,name))) 
  }
  
})

test_that("subset_scMethrix", {
  
  expect_error(subset_scMethrix("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_error(subset_scMethrix(scm),"At least 1 argument mandatory")
    
    samples <- colnames(scm)[c(1,3)] # Should be "C1" and "C3"
    contigs <- levels(seqnames(scm))[1] # Should be "chr1"
    regions <- GRanges(seqnames = levels(seqnames(scm))[1:2], ranges = IRanges(1,100000000)) 
    
    # Subset by include
    s <- subset_scMethrix(scm, samples = samples, by="include")
    expect_equal(dim(s),c(n_cpg,length(samples)))
    expect_equal(samples,colData(s)@rownames)
   
    s <- subset_scMethrix(scm, contigs = contigs, by="include")
    expect_equal(dim(s),c(length(which(as.vector(rowRanges(scm)@seqnames) %in% contigs)),n_samples))
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
    expect_equal(dim(s),c(n_cpg-length(which(as.vector(rowRanges(scm)@seqnames) %in% contigs)),n_samples))
    expect_equal(length(intersect(contigs,as.character(seqnames(s)@values))),0)
    
    s <- subset_scMethrix(scm, regions = regions, by = "exclude")
    expect_equal(dim(s),c(152,n_samples))
    expect_equal(length(findOverlaps(regions,rowRanges(s))),0)
    
    s <- subset_scMethrix(scm, samples = samples, contigs = contigs, regions = regions, by = "exclude")
    expect_equal(dim(s),c(72,length(samples)))
    
  }))
})

test_that("get_matrix", {
 
  expect_error(get_matrix("not scMethrix"),msg.validateExp)
  
    scm <- get_matrix(scm.h5)
    expect_equal(dim(scm),c(n_cpg,n_samples))  
    expect_is(scm,"HDF5Matrix")
    
    scm <- get_matrix(scm.mem)
    expect_equal(dim(scm),c(n_cpg,n_samples))  
    expect_is(scm,"matrix")
    
    invisible(lapply(list(scm.mem,scm.h5), function(scm) {
      expect_warning(get_matrix(scm,add_loci=FALSE, in_granges = TRUE))
      
      mtx <- get_matrix(scm=scm,add_loci=TRUE)
      expect_equal(dim(mtx),c(n_cpg,n_samples+3))  
      expect_is(mtx,"data.table")
      
      mtx <- get_matrix(scm,add_loci=TRUE, in_granges = TRUE)
     # expect_equal(seqnames(m)@lengths,c(10,8))
      expect_equal(dim(mcols(mtx)),c(n_cpg,n_samples))
      expect_is(mtx,"GRanges")
      
      mtx <- get_matrix(scm,order_by_sd = TRUE)
      expect_false(is.unsorted(rev(rowSds(mtx,na.rm=TRUE)),na.rm=TRUE))
      
      mtx <- get_matrix(scm,n_chunks = n_samples,by="col")
      invisible(lapply(mtx,function(m) expect_equal(dim(m),c(n_cpg,1))))
      
      mtx <- get_matrix(scm,n_chunks = n_cpg,by="row")
      invisible(lapply(mtx,function(m) expect_equal(dim(m),c(1,n_samples))))
      
  }))
})


test_that("remove_uncovered", {
  
  expect_error(remove_uncovered("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    samples = c("C1","C2")
    s <- subset_scMethrix(scm,samples=samples,by = "include")
    expect_equal(dim(s),c(n_cpg,length(samples)))
    non_na_rows = nrow(score(s)[rowSums(is.na(score(s))) != ncol(s), ])
    expect_equal(dim(remove_uncovered(s)),c(non_na_rows,length(samples)))
    
    expect_equal(remove_uncovered(s),remove_uncovered(s,n_threads=2))
    
  }))
})

test_that("get_stats", {
  
  expect_error(get_stats("not scMethrix"),msg.validateExp)
  
  fmt <- function(x) round(as.numeric(x),4)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    chrs <- seqlevels(scm)
    nChrs <- length(GenomeInfoDb::seqlengths(rowRanges(scm)))
    samples <- sampleNames(scm)
    nSamples <- nrow(colData(scm))
    mtx <- score(scm)
    
    # Check all
    stats <- get_stats(scm, per_sample = FALSE, per_chr = FALSE)
    expect_equal(dim(stats),c(1,6))
    expect_equal(fmt(mean(mtx, na.rm=TRUE)), fmt(stats$Mean))
    expect_equal(fmt(median(mtx,na.rm=TRUE)), fmt(stats$Median))
    expect_equal(fmt(sd(mtx,na.rm=TRUE)), fmt(stats$SD), tolerance=1e-1)
    expect_equal(sum(!is.na(mtx)), fmt(stats$Count))
    
    # Check per sample
    stats <- get_stats(scm, per_sample = TRUE, per_chr = FALSE)
    expect_equal(dim(stats),c(nSamples,6))
    expect_true(all(stats$Sample %in% samples))
    
    expect_equal(fmt(colMeans(mtx,na.rm=TRUE)), fmt(stats$Mean))
    expect_equal(fmt(colMedians(mtx,na.rm=TRUE)), fmt(stats$Median))
    expect_equal(fmt(colSds(mtx,na.rm=TRUE)), fmt(stats$SD), tolerance=1e-1)
    expect_equal(nrow(mtx)-fmt(colCounts(mtx,value=NA)), fmt(stats$Count))
    
    # Check per chr
    stats <- get_stats(scm, per_sample = FALSE, per_chr = TRUE)
    expect_equal(dim(stats),c(nChrs,6))
    expect_true(all(stats$Chromosome %in% chrs))
    
    chr <- stats$Chromosome[[1]]
    stats <- stats[1,]
    mtx <- score(subset_scMethrix(scm,contigs = chr))
    expect_equal(fmt(mean(mtx, na.rm=TRUE)), fmt(stats$Mean))
    expect_equal(fmt(median(mtx,na.rm=TRUE)), fmt(stats$Median))
    expect_equal(fmt(sd(mtx,na.rm=TRUE)), fmt(stats$SD), tolerance=1e-1)
    expect_equal(sum(!is.na(mtx)), fmt(stats$Count))
    
    # Check per sample, per chr
    stats <- get_stats(scm, per_sample = TRUE, per_chr = TRUE)
    expect_equal(dim(stats),c(nChrs*nSamples,6))
    expect_true(all(stats$Sample %in% samples))
    expect_true(all(stats$Chromosome %in% chrs))
    
    sample <- stats$Sample[[1]]
    chr <- stats$Chromosome[[1]]
    stats <- stats[1,]
    mtx <- score(subset_scMethrix(scm,contigs = chr, sample = sample))
    expect_equal(fmt(mean(mtx, na.rm=TRUE)), fmt(stats$Mean))
    expect_equal(fmt(median(mtx,na.rm=TRUE)), fmt(stats$Median))
    expect_equal(fmt(sd(mtx,na.rm=TRUE)), fmt(stats$SD), tolerance=1e-1)
    expect_equal(sum(!is.na(mtx)), fmt(stats$Count))
    
    # Check for stat subset
    stats <- get_stats(scm,stat="Mean")
    expect_equal(dim(stats),c(nChrs*nSamples,3))
    expect_true(colnames(stats)[ncol(stats)] == "Mean")
    
    stats <- get_stats(scm,stat="Median")
    expect_equal(dim(stats),c(nChrs*nSamples,3))
    expect_true(colnames(stats)[ncol(stats)] == "Median")
    
    stats <- get_stats(scm,stat="Count")
    expect_equal(dim(stats),c(nChrs*nSamples,3))
    expect_true(colnames(stats)[ncol(stats)] == "Count")
    
    stats <- get_stats(scm,stat="SD")
    expect_equal(dim(stats),c(nChrs*nSamples,3))
    expect_true(colnames(stats)[ncol(stats)] == "SD")
    
    stats <- get_stats(scm,stat=c("Mean","Median"))
    expect_equal(dim(stats),c(nChrs*nSamples,4))
    expect_true(any(colnames(stats) == "Mean"))
    expect_true(any(colnames(stats) == "Median"))
    
    # Check ignores
    stats <- get_stats(scm, per_chr = FALSE, ignore_sample = samples[-(1)])
    expect_equal(stats$Sample,samples[1])
    
    stats <- get_stats(scm, per_sample = FALSE, ignore_chr = chrs[-(1)])
    expect_equal(stats$Chromosome,chrs[1])
  }))
})

test_that("get_region_summary", {
  
  expect_error(get_region_summary("not scMethrix"),msg.validateExp)
  
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
# 
# test_that("mask_by_coverage", {
#   expect_error(mask_by_coverage("not scMethrix"),msg.validateExp)
#   expect_error(mask_by_coverage(scm.mem,assay="not an assay"),msg.validateAssay)
#   expect_error(mask_by_coverage(remove_assay(scm.mem,assay="counts")))
#   expect_error(mask_by_coverage(scm.mem,n_threads=2))
#   expect_error(mask_by_coverage(scm.mem,low_threshold=-1),"low_threshold")
#   expect_error(mask_by_coverage(scm.mem,low_threshold="not numeric"),msg.validateType)
#   expect_error(mask_by_coverage(scm.mem,avg_threshold=-1,"avg_threshold"))
#   expect_error(mask_by_coverage(scm.mem,avg_threshold="not numeric"),msg.validateType)
#   
#   invisible(lapply(list(scm.mem,scm.h5), function(scm) {
# 
#     scm <- mask_by_coverage(scm,low_threshold=2,avg_threshold=NULL)
#     expect_equal(dim(scm),c(n_cpg,n_samples))
#     expect_equal(dim(remove_uncovered(scm)),c(length(which(rowSums(counts(scm),na.rm=TRUE) >= 2)),n_samples))
# 
#     scm <- mask_by_coverage(scm,low_threshold=NULL,avg_threshold=1) #Removes all rows with a sample with a coverage of 2
#     expect_equal(dim(scm),c(n_cpg,n_samples))
#     expect_equal(dim(remove_uncovered(scm)),c(length(which(rowMeans(counts(scm),na.rm=TRUE) == 1)),n_samples))
# 
#   }))
# })
# 
# test_that("mask_by_sample", {
#   expect_error(mask_by_sample("not scMethrix"),msg.validateExp)
#   expect_error(mask_by_sample(scm.mem,assay="not an assay"),msg.validateAssay)
#   expect_error(mask_by_sample(scm.mem,n_threads=2))
#   expect_error(mask_by_sample(scm.mem,low_threshold=2,prop_threshold=1))
#   expect_error(mask_by_sample(scm.mem,low_threshold=-1),"low_threshold")
#   expect_error(mask_by_sample(scm.mem,low_threshold="not numeric"),msg.validateType)
#   expect_error(mask_by_sample(scm.mem,prop_threshold=-1,"prop_threshold"))
#   expect_error(mask_by_sample(scm.mem,prop_threshold=2,"prop_threshold"))
#   expect_error(mask_by_sample(scm.mem,prop_threshold="not numeric"),msg.validateType)
#   
#   invisible(lapply(list(scm.mem,scm.h5), function(scm) {
#     
#     scm <- mask_by_sample(scm,low_threshold=2)
#     expect_equal(dim(scm),c(n_cpg,n_samples))
#     expect_equal(dim(remove_uncovered(scm)),c(length(which(ncol(scm) - rowCounts(score(scm),value=NA) >= 2)),n_samples))
#     
#     scm <- mask_by_sample(scm,low_threshold=NULL,prop_threshold=0.25) # Since 0.25 of 4 sample is 1, same result as low_threshold = 1
#     expect_equal(dim(scm),c(n_cpg,n_samples))
#     expect_equal(dim(remove_uncovered(scm)),c(length(which(ncol(scm) - rowCounts(score(scm),value=NA) > 1)),n_samples))
#     
#   }))
# })
# 
# test_that("mask_by_variance", {
#   expect_error(mask_by_variance("not scMethrix"),msg.validateExp)
#   expect_error(mask_by_variance(scm.mem,assay="not an assay"),msg.validateAssay)
#   expect_error(mask_by_variance(scm.mem,n_threads=2))
#   expect_error(mask_by_variance(scm.mem,low_threshold=2,"low_threshold must be between 0 and 1"))
#   expect_error(mask_by_variance(scm.mem,low_threshold=-1,"low_threshold must be between 0 and 1"))
#   expect_error(mask_by_variance(scm.mem,low_threshold="not numeric"),msg.validateType)
#   
#   invisible(lapply(list(scm.mem,scm.h5), function(scm) {
#     
#     #Since only 4 samples, only rows with 1 unique score value are masked
#     uniq <- apply(score(scm),1,function(x) {
#       uniq <- unique(x)
#       return(length(uniq[!is.na(uniq)]))
#     })
#     
#     scm <- mask_by_variance(scm,low_threshold=0.05)
# 
#     expect_equal(dim(scm),c(n_cpg,n_samples))
#     expect_equal(dim(remove_uncovered(scm)),c(sum(uniq != 1),n_samples))
#   }))
# })

test_that("mask_scMethrix", {

  expect_error(mask_scMethrix("not scMethrix"),msg.validateExp)
  expect_error(mask_scMethrix(scm.mem,n_threads=2))

  invisible(lapply(list(scm.mem,scm.h5), function(scm) {

    expect_error(mask_scMethrix(scm,assay="not an assay"),msg.validateAssay)
    expect_error(mask_scMethrix(scm,by="not an arg"),msg.validateArg)
    expect_error(mask_scMethrix(scm,stat="not an arg"),msg.validateArg)
    expect_error(mask_scMethrix(scm,op="not an arg"),msg.validateArg)
    expect_error(mask_scMethrix(scm,na.rm="not a boolean"),msg.validateType)
    expect_error(mask_scMethrix(scm,verbose="not a boolean"),msg.validateType)
    expect_error(mask_scMethrix(scm,threshold="not a numeric"),msg.validateType)

    msg.cpgErr = "No CpG sites left"

    s <- mask_scMethrix(scm,assay="counts",threshold=2,by="row",stat="sum",op="==")
    row_idx <- which(DelayedMatrixStats::rowSums2(counts(scm),na.rm=T) == 2)
    expect_true(all(is.na(score(s[row_idx,]))))

    s <- mask_scMethrix(scm,assay="counts",threshold=1,by="row",stat="mean",op=">")
    row_idx <- which(DelayedMatrixStats::rowMeans2(counts(scm)) > 1)
    expect_true(all(is.na(score(s[row_idx,]))))

    avg <- nrow(scm) - mean(DelayedMatrixStats::colCounts(score(scm),na.rm=T,value = as.integer(NA)))
    s <- mask_scMethrix(scm,assay="score",threshold=avg,by="col",stat="count",op=">")
    col_idx <- (nrow(scm) - DelayedMatrixStats::colCounts(score(scm),na.rm=T,value = as.integer(NA))) > avg
    expect_true(all(is.na(score(s[,col_idx]))))

    s <- remove_uncovered(scm)
    s <- mask_scMethrix(s,assay="counts",threshold=0,by="row",stat="sum",op="==")
    expect_false(any(DelayedMatrixStats::rowAlls(score(s),value=NA)))
    
    s <- remove_uncovered(scm)
    expect_error(mask_scMethrix(s,assay="counts",threshold=0,by="row",stat="sum",op=">"),msg.cpgErr)
    
    
    expect_error(mask_scMethrix(scm,assay="score",threshold=0,by="row",stat="var",op=">="),msg.cpgErr)
    expect_error(mask_scMethrix(scm,assay="score",threshold=1,by="row",stat="var",op="<="),msg.cpgErr)

  }))
})

