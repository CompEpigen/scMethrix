test_that("read_index", {

  index <- read_index(files, col_list = col_list)
  expect_equal(names(index),c("chr","start","end"))
  expect_equal(dim(index),c(n_cpg,3))
  expect_equal(index,read_index(files,n_threads=2, col_list = col_list))
  
})

test_that("read_bed_by_index", {
  
  file <- files[1]
  index <- read_index(files, col_list = col_list)
  C1 <- data.table::fread(file = file, sep = '\t', header = FALSE)
  bed <- read_bed_by_index(file,index,col_list=col_list)[[1]]
  expect_equal(dim(bed),c(nrow(index),1))
  expect_equal(dim(bed),c(nrow(C1),1))

})

test_that("read_bed - Input errors", {
  
  expect_error(read_beds(NULL))
  expect_error(read_beds(tempfile))
  
})

test_that("read_bed - HDF5, no coverage", {
  
  #expect_error(read_beds(files,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, h5=TRUE,h5_dir=NULL))
  
  path <- paste0(h5_dir,"/HDF5mem")
  unlink(path, recursive = TRUE)
  suppressWarnings(dir.create(path,recursive=TRUE))
  
  scm1 <- read_beds(files, is_h5=TRUE, h5_dir=path, replace=TRUE, chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, colData = colData)
  scm2 <- load_scMethrix(dest=path)
    
  expect_true(is_h5(scm1))
  #expect_false(has_cov(scm1))
  expect_equal(class(scm1)[1], "scMethrix")
  expect_identical(names(assays(scm1)), "score")
  expect_identical(rownames(colData(scm1)),rownames(colData))
  expect_equal(class(get_matrix(scm1))[[1]], "HDF5Matrix")
  expect_equal(dim(scm1),c(n_cpg, n_samples))
  expect_equivalent(scm1, scm2)
  
  unlink(path, recursive = TRUE)
  
})

test_that("read_bed - HDF5, with coverage", {
  
  # expect_error(read_beds(files,h5=TRUE,h5_dir=NULL,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5),
  #              msg.validateType)
  # 
  path <- paste0(h5_dir,"/HDF5mem")
  unlink(path, recursive = TRUE)
  suppressWarnings(dir.create(path,recursive=TRUE))
  
  scm1 <- read_beds(files, is_h5=TRUE, h5_dir=path, replace=TRUE, chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)
  scm2 <- load_scMethrix(dest=path)
  
  expect_true(is_h5(scm1))
  expect_true(has_cov(scm1))
  expect_equal(class(scm1)[1],"scMethrix")
  expect_identical(names(assays(scm1)),c("score","counts"))
  expect_identical(rownames(colData(scm1)),rownames(colData))
  expect_equal(class(get_matrix(scm1))[[1]],"HDF5Matrix")
  expect_equal(class(get_matrix(scm1,assay="counts"))[[1]],"HDF5Matrix")
  expect_equal(dim(scm1),c(n_cpg, n_samples))
  expect_equivalent(scm1,scm2)
  
  unlink(path, recursive = TRUE)
  
})

test_that("read_bed - in-memory, no coverage", {
  
  scm <- read_beds(files,is_h5=FALSE, chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, colData = colData)
  
  expect_equal(class(scm)[1],"scMethrix")
  expect_identical(names(assays(scm)),"score")
  expect_identical(rownames(colData(scm)),rownames(colData))
  expect_equal(class(get_matrix(scm))[[1]],"matrix")
  expect_equal(dim(scm),c(n_cpg, n_samples))
  
})

test_that("read_bed - in-memory, with coverage", {
  
  scm <- read_beds(files,is_h5=FALSE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)
  
  expect_equal(class(scm)[1],"scMethrix")
  expect_identical(names(assays(scm)),c("score","counts"))
  expect_identical(rownames(colData(scm)),rownames(colData))
  expect_equal(class(get_matrix(scm))[[1]],"matrix")
  expect_equal(class(get_matrix(scm,assay="counts"))[[1]],"matrix")
  expect_equal(dim(scm),c(n_cpg, n_samples))
  
})

test_that("read_bed - HDF5 and in-memory equivalence, no coverage", {
  
  path <- paste0(h5_dir,"HDF5equiv")
  
  h5 <- read_beds(files,is_h5=TRUE,h5_dir=path,replace=TRUE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, colData = colData)
  mem <- read_beds(files,is_h5=FALSE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, colData = colData)
  
  expect_true(identical.scm(h5,mem,exclude_is_h5=TRUE))
  
  # expect_equal(as.matrix(score(scm.hdf)),score(scm.mem))
  # expect_equal(rowRanges(scm.hdf),rowRanges(scm.mem))
  
  unlink(path, recursive = TRUE)
})

test_that("read_bed - HDF5 and in-memory equivalence, with coverage", {
  
  path <- paste0(h5_dir,"HDF5equiv")
  
  h5 <- read_beds(files,is_h5=TRUE,h5_dir=path,replace=TRUE, chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)
  mem <- read_beds(files,is_h5=FALSE, chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)

  expect_true(identical.scm(h5,mem,exclude_is_h5=TRUE))
  
  # expect_equal(as.matrix(score(scm.hdf)),score(scm.mem))
  # expect_equal(as.matrix(counts(scm.hdf)),counts(scm.mem))
  # expect_equal(rowRanges(scm.hdf),rowRanges(scm.mem))
  
  unlink(path, recursive = TRUE)
})

test_that("read_bed - threaded", {

  path <- paste0(h5_dir,"threaded")
  
    scm <- lapply(c(0,2), function(x) read_beds(files,is_h5=FALSE,n_threads = x, chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5))
    expect_true(validObject(scm))
    expect_equal(dim(scm[[1]]),c(n_cpg,4))
    expect_equal(scm[[1]],scm[[2]])
    
    scm1 <- read_beds(files,is_h5=TRUE,h5_dir=paste0(h5_dir,"1"),replace=TRUE,n_threads = 1,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5)
    scm2 <- read_beds(files,is_h5=TRUE,h5_dir=paste0(h5_dir,"2"),replace=TRUE,n_threads = 2,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5)
    expect_equal(dim(scm[[1]]),c(n_cpg,4))
    expect_equal(as.matrix(score(scm[[1]])),as.matrix(score(scm[[2]])))
    expect_equal(as.matrix(counts(scm[[1]])),as.matrix(counts(scm[[2]])))
    expect_equal(metadata(scm[[1]]),metadata(scm[[2]]))
  
})

test_that("read_bed - batched", {
  
  path <- paste0(h5_dir,"batched")
  
  scm <- read_beds(files,is_h5=TRUE, h5_dir=path,replace=TRUE, chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, batch_size = n_samples-1)

  expect_equal(as.matrix(score(scm)),as.matrix(score(scm.h5)))
  expect_equal(as.matrix(counts(scm)),as.matrix(counts(scm.h5)))
  
  scm <- read_beds(files,is_h5=FALSE, chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, batch_size = n_samples-1)
  
  expect_equal(as.matrix(score(scm)),as.matrix(score(scm.mem)))
  expect_equal(as.matrix(counts(scm)),as.matrix(counts(scm.mem)))

})
