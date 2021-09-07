test_that("get_distance_matrix", {
  
  types = c("pearson", "spearman", "tau", "euclidean", "maximum", 
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
    
    expect_error(get_distance_matrix(scm,assay="impute",type="not a metric"),"Invalid type of distance calculation")
    
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
  
  types = c("hierarchical", "partition", "model")
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
    
    expect_error(cluster_scMethrix(scm,assay="impute",type="not a type"), "Invalid type of clustering")
    
    dist <- get_distance_matrix(scm, assay="impute")

    expect_equivalent(ncol(colData(scm)),0) # Check there's no colData before clustering
    
    invisible(lapply(types, function(type) {

      # Looks weird, but should throw warning if running the same cluster_scMethrix command back to back (replaces the colData "Cluster" column)
      if (type == "model") {
        expect_warning(
          scm <- cluster_scMethrix(scm,n_clusters=n_clusters,assay="impute",type=type,colname=name)
        ,"n_clusters is ignored")
        expect_warning(
          scm <- cluster_scMethrix(scm,n_clusters=n_clusters,assay="impute",type=type,colname=name)
          ,"Colnames of colData already exist")
      } else {
        scm <- cluster_scMethrix(scm,n_clusters=n_clusters,assay="impute",type=type,colname=name)
        expect_warning(scm <- cluster_scMethrix(scm,n_clusters=n_clusters,assay="impute",type=type,colname=name),"Colnames of colData already exist")
      }
      
      cd <- colData(scm)
      expect_is(scm,"scMethrix")
      expect_equivalent(colnames(cd),name)
      expect_true(all(cd$Cluster %in% 1:n_clusters))
    }))
  }))
})

test_that("append_col_data", {

  expect_error(append_col_data(scm="not scMethrix"),"A valid scMethrix object needs to be supplied")
  #expect_error(append_col_data(scm=scm.mem, colData=NULL),"A valid colData object must be supplied")
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_equivalent(ncol(colData(scm)),0) # Check there's no colData before appending
    
    name = "Cluster"
    vals <- 1:ncol(scm)
    
    #Dataframe input
    colData <- colData(scm)
    colData[name] <- vals
    app <- append_col_data(scm, colData=colData,name=name)
    
    expect_equivalent(colnames(colData(app)),name)
    expect_equivalent(colData(app)[,name],vals)
    
    #Named vector input
    colData <- vals
    names(colData) <- rownames(colData(scm))
    app <- append_col_data(scm, colData=colData,name=name)
    
    expect_equivalent(colnames(colData(app)),name)
    expect_equivalent(colData(app)[,name],vals)
  }))
})