test_that("read_index", {

  index <- read_index(files)
  expect_equivalent(names(index),c("chr","start","end"))
  expect_equivalent(dim(index),c(n_cpg,3))
  expect_equivalent(read_index(files),read_index(files,n_threads=2))
  
})

test_that("read_bed_by_index", {
  
  file <- files[1]
  index <- read_index(files)
  C1 <- read.table(file = file, sep = '\t', header = FALSE)
  bed <- read_bed_by_index(file,index)[[1]]
  expect_equivalent(dim(bed),c(nrow(index),1))
  expect_equivalent(dim(na.omit(bed)),c(nrow(C1),1))

})

test_that("read_bed - Input errors", {
  
  expect_error(read_beds(NULL))
  expect_error(read_beds(tempfile))
  
})

test_that("read_bed - HDF5, no coverage", {
  
  expect_error(read_beds(files,h5=TRUE,h5_dir=NULL))
  
  path <- paste0(h5_dir,"HDF5mem")
  
  scm1 <- read_beds(files,h5=TRUE,h5_dir=path,replace=TRUE)
  scm2 <- load_HDF5_scMethrix(dir=path)
    
  expect_true(is_h5(scm1))
  expect_false(has_cov(scm1))
  expect_equivalent(class(scm1)[1],"scMethrix")
  expect_equivalent(class(get_matrix(scm1))[[1]],"HDF5Matrix")
  expect_equivalent(dim(scm1),c(n_cpg,4))
  expect_equivalent(scm1,scm2)
  
  unlink(path, recursive = TRUE)
  
})

test_that("read_bed - in-memory, no coverage", {
  
  scm <- read_beds(files,h5=FALSE)
  
  expect_equivalent(class(scm)[1],"scMethrix")
  expect_equivalent(class(get_matrix(scm))[[1]],"matrix")
  expect_equivalent(dim(scm),c(n_cpg,4))
  
})

test_that("read_bed - HDF5, with coverage", {
  
  expect_error(read_beds(files,h5=TRUE,h5_dir=NULL))
  
  path <- paste0(h5_dir,"HDF5mem")
  
  scm1 <- read_beds(files,h5=TRUE,h5_dir=path,replace=TRUE,cov_idx=5)
  scm2 <- load_HDF5_scMethrix(dir=path)
  
  expect_true(is_h5(scm1))
  expect_true(has_cov(scm1))
  expect_equivalent(class(scm1)[1],"scMethrix")
  expect_equivalent(class(get_matrix(scm1))[[1]],"HDF5Matrix")
  expect_equivalent(class(get_matrix(scm1,assay="counts"))[[1]],"HDF5Matrix")
  expect_equivalent(dim(scm1),c(n_cpg,4))
  expect_equivalent(scm1,scm2)
  
  unlink(path, recursive = TRUE)
  
})

test_that("read_bed - in-memory, with coverage", {
  
  scm <- read_beds(files,h5=FALSE,cov_idx=5)
  
  expect_equivalent(class(scm)[1],"scMethrix")
  expect_equivalent(class(get_matrix(scm))[[1]],"matrix")
  expect_equivalent(class(get_matrix(scm,assay="counts"))[[1]],"matrix")
  expect_equivalent(dim(scm),c(n_cpg,4))
  
})

test_that("read_bed - HDF5 and in-memory equivalence, no coverage", {
  
  path <- paste0(h5_dir,"HDF5equiv")
  
  scm.hdf <- read_beds(files,h5=TRUE,h5_dir=path,replace=TRUE)
  scm.mem <- read_beds(files,h5=FALSE)
  
  expect_equivalent(as.matrix(score(scm.hdf)),score(scm.mem))
  expect_equivalent(rowRanges(scm.hdf),rowRanges(scm.mem))
  
  unlink(path, recursive = TRUE)
})

test_that("read_bed - HDF5 and in-memory equivalence, with coverage", {
  
  path <- paste0(h5_dir,"HDF5equiv")
  
  scm.hdf <- read_beds(files,h5=TRUE,h5_dir=path,replace=TRUE, cov_idx = 5)
  scm.mem <- read_beds(files,h5=FALSE, cov_idx = 5)
  
  expect_equivalent(as.matrix(score(scm.hdf)),score(scm.mem))
  expect_equivalent(as.matrix(counts(scm.hdf)),counts(scm.mem))
  expect_equivalent(rowRanges(scm.hdf),rowRanges(scm.mem))
  
  unlink(path, recursive = TRUE)
})

test_that("read_bed - threaded", {
  
  scm <- lapply(c(0,2), function(x) read_beds(files,h5=FALSE,n_threads = x))
  expect_equivalent(dim(scm[[1]]),c(n_cpg,4))
  expect_equivalent(scm[[1]],scm[[2]])
  
  scm <- lapply(c(0,2), function(x) read_beds(files,h5=FALSE,cov_idx=5,n_threads = x))
  expect_equivalent(dim(scm[[1]]),c(n_cpg,4))
  expect_equivalent(scm[[1]],scm[[2]])
  
  # scm <- lapply(c(0,2), function(x) read_beds(files,h5=TRUE,h5_dir=h5_dir,replace=TRUE,n_threads = x))
  # expect_equivalent(dim(scm[[1]]),c(100,4))
  # expect_equivalent(scm[[1]],scm[[2]])
  # 
  # scm <- lapply(c(0,2), function(x) read_beds(files,h5=TRUE,h5_dir=h5_dir,replace=TRUE,cov_idx=5,n_threads = x))
  # expect_equivalent(dim(scm[[1]]),c(100,4))
  # expect_equivalent(scm[[1]],scm[[2]])

})
