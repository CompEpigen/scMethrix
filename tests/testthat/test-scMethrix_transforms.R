test_that("bin_scMethrix", {

  expect_error(bin_scMethrix("not scMethrix"),msg.validateExp)
  expect_error(bin_scMethrix(scm.h5,regions=GRanges()),msg.validateType)

  path <- paste0(h5_dir,"bin")
  regions <- GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(1,1000000000)) 
  
 invisible(lapply(list(scm.mem,scm.h5), function(scm) {
   
    # Check default conditions and threading
    bin <- bin_scMethrix(scm, h5_dir = paste0(h5_dir,"/bin1"), n_threads = 2, replace = T)
    expect_equal(dim(bin),c(256,4))
    
    #Check the score assay (should be mean)
    scm <- transform_assay(scm,assay="score",new_assay="bin",trans=binarize)
    scm <- transform_assay(scm,assay="score",new_assay="bin2",trans=binarize)
    bin <- bin_scMethrix(scm,bin_size=1000,bin_by="cpg", h5_dir = paste0(h5_dir,"/bin1"), replace = T)
    expect_equal(dim(bin),c(length(rowRanges(bin)),length(colnames(bin))))
    
    if (is_h5(scm)) {
      expect_equal(class(score(bin))[[1]],"DelayedMatrix")
    }
    
    sub <- subset_scMethrix(scm,contigs="chr1")
    vals <- DelayedMatrixStats::colMeans2(score(sub),na.rm=T)
    expect_equal(as.numeric(score(bin)[1,]),as.numeric(vals))
    
    #Check the counts assay (should be sum)
    vals <- DelayedMatrixStats::colSums2(counts(sub),na.rm=T)
    expect_equal(as.numeric(counts(bin)[1,]),as.numeric(vals))
    
    #Check a custom assay  (should be mean)
    vals <- DelayedMatrixStats::colMeans2(get_matrix(sub,assay="bin"),na.rm=T)
    expect_equal(as.numeric(get_matrix(bin,assay="bin")[1,]),as.numeric(vals))
    
    #Check the custom transform function  (should be mean, but specified as sum)
    expect_error(bin_scMethrix(scm,trans="Not a trans"))
    
    bin2 <- bin_scMethrix(scm,bin_size=1000,bin_by="cpg",trans = c(bin2 = function(x) sum(x,na.rm=TRUE)), 
                          h5_dir = paste0(h5_dir,"/bin2"),replace=T)
    expect_equal(score(bin),score(bin2), check.attributes = FALSE)
    expect_equal(counts(bin),counts(bin2), check.attributes = FALSE)
    expect_equal(get_matrix(bin,assay="bin"),get_matrix(bin2,assay="bin"), check.attributes = FALSE)
    
    vals <- DelayedMatrixStats::colSums2(get_matrix(sub,assay="bin2"),na.rm=T)
    expect_equal(as.numeric(get_matrix(bin2,assay="bin2")[1,]),as.numeric(vals))
    
    #Check for region subsetting
    bin <- bin_scMethrix(scm,regions=regions,bin_size=1000,bin_by="cpg", h5_dir = paste0(h5_dir,"/bin3"), replace = T)
    expect_equal(sum(rowRanges(bin)$n_cpgs),n_cpg)
  }))
})

test_that("collapse_samples", {
  
  expect_error(collapse_samples("not scMethrix"),msg.validateExp)
  
  #invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  invisible(lapply(list(scm.mem), function(scm) {
  
    groups <- sort(rep_len(c("Grp1","Grp2"),nrow(colData(scm))))
    
    colData(scm)["Cluster"] <- groups #TODO Should be dynamic in case input data changes
    
    expect_error(collapse_samples(scm, colname="Not a colname"),"Cannot find column")
    expect_error(collapse_samples(scm, colname = "Cluster",trans="not a trans"))
    
    expect_is(scm.col <- collapse_samples(scm, colname = "Cluster"),"scMethrix")
    
    expect_equal(sort(sampleNames(scm.col)),sort(unique(colData(scm)$Cluster)))
    
    grp = groups[1]
    grp_len = sum(groups == grp)
    
    expect_equal(as.vector(score(scm.col)[,grp]),rowMeans(score(scm)[,1:grp_len],na.rm=TRUE))
    expect_equal(as.vector(counts(scm.col[,grp])),rowSums(counts(scm)[,1:grp_len],na.rm=TRUE))
    expect_equal(colData(scm.col)$n_Samples,as.vector(table(groups)))
    
  }))
})

test_that("transform_assay", {
  
  expect_error(transform_assay("not scMethrix"),msg.validateExp)
  trans <- function(x) x+1
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(transform_assay(scm,trans="not closure"),msg.validateType)
    expect_warning(transform_assay(scm, trans=trans, assay="score",new_assay="score"))
    
    # Create a new assay with value of x+1
    plus1 <- transform_assay(scm, trans=trans, assay="score",new_assay="plus1")
    expect_false(isTRUE(all.equal(assays(scm), assays(plus1))))
    expect_equivalent(score(scm)+1,get_matrix(plus1,assay="plus1"))
    if (is_h5(scm)) {
      expect_equal(class(get_matrix(plus1,assay="plus1"))[[1]],"HDF5Matrix")
    }
    
    #Make sure the other assays are not affected
    expect_equivalent(counts(scm),counts(plus1))
  }))
})

test_that("impute_regions", {
  expect_error(impute_regions("not scMethrix"),msg.validateExp)
  expect_warning(impute_regions(scm.h5),"Imputation cannot be done on HDF5 data. Data will be cast as matrix for imputation.")

  fun <- function(mtx) missForest::missForest(mtx)$ximp
  
  suppressWarnings(
    # Check all the usable imputation methods
    lapply(list("kNN","iPCA","RF",fun), function (method) {
    invisible(lapply(list(scm.mem,scm.h5), function(scm) {
        expect_error(impute_regions(scm,assay = "not an assay"),msg.validateAssay)
        expect_error(impute_regions(scm,new_assay = "score"))
        expect_warning(impute_regions(scm,new_assay = "counts",type=method))
        
        impute = impute_regions(scm,new_assay="impute",type=method)
        expect_true("impute" %in% SummarizedExperiment::assayNames(impute))
        
        sco <- get_matrix(impute,assay="score")
        imp <- get_matrix(impute,assay="impute")
        NAs <- which(is.na(sco))
        nonNAs <- which(!is.na(sco))
        
        expect_true(anyNA(sco) && !anyNA(imp))
        expect_equivalent(sco[nonNAs],imp[nonNAs])
        expect_false(all(sco[NAs] %in% imp[NAs]))
      }))
    })
  )
  
  expect_error(impute_regions(scm.mem, type=function(x) sum(x)),"Error with imputation algorithm")
  
})

test_that("generate_training_set", {
  
  expect_error(generate_training_set("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(generate_training_set(scm,training_prop = 2),"training_prop must in the range of")
    
    n <- nrow(scm)
    prop <- 0.2
    
    set <- generate_training_set(scm,training_prop = prop)
    expect_equal(nrow(set$training),round(n*prop))
    expect_equal(nrow(set$test),round(n*(1-prop)))
    expect_equivalent(scm,merge_scMethrix(set$training,set$test))
  }))
})

test_that("generate_random_subset", {
  
  expect_error(generate_random_subset("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_warning(generate_random_subset(scm,n_cpgs = nrow(scm)+1))
    expect_warning(generate_random_subset(scm,n_cpgs = -1))
    
    m <- generate_random_subset(scm,n_cpgs = 100)
    expect_equal(nrow(m),100)
  }))
})
