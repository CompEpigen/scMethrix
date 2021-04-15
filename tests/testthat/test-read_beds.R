test_that("read_index", {

  expect_equivalent(read_index(files),rbind(df1,df3,df4)[,1:3])
  #expect_equivalent(read_index(files),read_index(files,n_threads=2))
  
})

test_that("read_bed_by_index", {
  
  index <- read_index(files)
  expect_equivalent(as.vector(unlist(read_bed_by_index(files[1],index)[1])),c(rep(0,5),rep(NA,13)))
  
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
  expect_equivalent(dim(scm1),c(18,4))
  expect_equivalent(scm1,scm2)
  
  unlink(path, recursive = TRUE)
  
})

test_that("read_bed - in-memory, no coverage", {
  
  scm <- read_beds(files,h5=FALSE)
  
  expect_equivalent(class(scm)[1],"scMethrix")
  expect_equivalent(class(get_matrix(scm))[[1]],"matrix")
  expect_equivalent(dim(scm),c(18,4))
  
})

test_that("read_bed - HDF5 and in-memory equivalence, no coverage", {

  path <- paste0(h5_dir,"HDF5equiv")
  
  scm.hdf <- read_beds(files,h5=TRUE,h5_dir=path,replace=TRUE)
  scm.mem <- read_beds(files,h5=FALSE)
  
  expect_equivalent(as.matrix(assays(scm.hdf)$score),assays(scm.mem)$score)
  expect_equivalent(rowRanges(scm.hdf),rowRanges(scm.mem))
  
  unlink(path, recursive = TRUE)
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
  expect_equivalent(class(get_matrix(scm1,type="c"))[[1]],"HDF5Matrix")
  expect_equivalent(dim(scm1),c(18,4))
  expect_equivalent(scm1,scm2)
  
  unlink(path, recursive = TRUE)
  
})
