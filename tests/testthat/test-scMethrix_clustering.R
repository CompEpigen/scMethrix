test_that("get_distance_matrix", {
  
  types = c("pearson", "spearman", "kendall", "euclidean", 
            "manhattan", "canberra", "binary", "minkowski")
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_error(get_distance_matrix(scm="not scMethrix"),"A valid scMethrix object needs to be supplied")
    expect_error(get_distance_matrix(scm),"There are NA values present")
    expect_error(get_distance_matrix(scm,assay="not an assay"),"Assay does not exist")
    
    if (is_h5(scm)) {
      expect_warning(scm <- impute_regions(scm), "Imputation cannot be done on HDF5 data.")
    } else {
      scm <- impute_regions(scm) 
    }
    
    expect_error(get_distance_matrix(scm,assay="impute",type="not a metric"),"Invalid input")
    
    invisible(lapply(types, function(metric) {
      dist <- get_distance_matrix(scm, assay="impute",type=metric)
      expect_equal(dim(as.matrix(dist)),rep(ncol(scm),2))
      expect_equal(row.names(as.matrix(dist)),row.names(colData(scm)))
      expect_false(any(is.na(as.matrix(dist))))
      return(dist)
    }))
  }))
})

test_that("cluster_scMethrix", {
  
  types = c("hier", "part", "model")
  name = "Cluster"
  n_clusters = 3
  
  expect_error(cluster_scMethrix(scm="not scMethrix"),"A valid scMethrix object needs to be supplied")
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_error(cluster_scMethrix(scm="not scMethrix"),"A valid scMethrix object needs to be supplied")
    expect_error(cluster_scMethrix(scm,assay="not an assay"),"Assay does not exist")
    
    if (is_h5(scm)) {
      expect_warning(scm <- impute_regions(scm), "Imputation cannot be done on HDF5 data.")
    } else {
      scm <- impute_regions(scm) 
    }
    
    expect_error(cluster_scMethrix(scm,assay="impute",type="not a type"), "Invalid input")
    
    #dist <- get_distance_matrix(scm, assay="impute")

    #expect_equivalent(ncol(colData(scm)),0) # Check there's no colData before clustering
    
    invisible(lapply(types, function(type) {
      scm.c <- scm
      if (type == "model") {
        expect_warning(
          scm.c <- cluster_scMethrix(scm.c,n_clusters=n_clusters,assay="impute",type=type,colname=name)
          ,"n_clusters is ignored")
      } else {
        scm.c <- cluster_scMethrix(scm.c,n_clusters=n_clusters,assay="impute",type=type,colname=name)
      }
      
      cd <- colData(scm.c)
      expect_is(scm.c,"scMethrix")
      expect_equivalent(colnames(cd),c(colnames(colData(scm)),name))
      expect_true(all(cd$Cluster %in% 1:n_clusters))
    }))
  }))
})

test_that("append_colData", {

  expect_error(append_colData(scm="not scMethrix"),"A valid scMethrix object needs to be supplied")
  #expect_error(append_colData(scm=scm.mem, colData=NULL),"A valid colData object must be supplied")
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
   # expect_equivalent(ncol(colData(scm)),0) # Check there's no colData before appending

    name = "Cluster"
    vals <- 1:ncol(scm)
    
    #Dataframe input
    colData <- colData(scm)
    colData[name] <- vals
    colData <- subset(colData, select=name)
    app <- append_colData(scm, colData=colData)
    expect_warning(app <- append_colData(app, colData=colData),"Colnames of colData already exist")
    expect_equivalent(colnames(colData(app)),c(colnames(colData(scm)),name))
    expect_equivalent(colData(app)[,name],vals)

    #Named vector input
    colData <- vals
    names(colData) <- rownames(colData(scm))
    app <- append_colData(scm, colData=colData,name=name)
    expect_warning(app <- append_colData(app, colData=colData,name=name),"Colnames of colData already exist")
    expect_equivalent(colnames(colData(app)),c(colnames(colData(scm)),name))
    expect_equivalent(colData(app)[,name],vals)

  }))
})