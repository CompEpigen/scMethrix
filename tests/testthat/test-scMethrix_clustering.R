test_that("get_distance_matrix", {
  
  types = c("pearson", "spearman", "kendall", "euclidean", 
            "manhattan", "canberra", "binary", "minkowski")
  
  expect_error(get_distance_matrix("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_error(get_distance_matrix(scm="not scMethrix"),msg.validateExp)
    expect_error(get_distance_matrix(scm),"There are NA values present")
    expect_error(get_distance_matrix(scm,assay="not an assay"),msg.validateAssay)
    
    scm <- transform_assay(scm,new_assay="fill", trans = function(x) fill(x,fill=0))
    
    expect_error(get_distance_matrix(scm,assay="fill",type="not a metric"),msg.validateArg)
    
    invisible(lapply(types, function(metric) {
      
      if (is_h5(scm)) {
        expect_warning(dist <- get_distance_matrix(
          scm, assay="fill",type=metric), "Distance matrix cannot be generated for HDF5 data")
      } else {
        dist <- get_distance_matrix(scm, assay="fill",type=metric)
      }
      
      expect_equal(dim(as.matrix(dist)),rep(ncol(scm),2))
      expect_equal(row.names(as.matrix(dist)),sampleNames(scm))
      expect_false(any(is.na(as.matrix(dist))))
      return(dist)
    }))
  }))
})

test_that("cluster_scMethrix", {
  
  types = c("hier", "part", "model")
  name = "Cluster"
  n_clusters = 3
  
  expect_error(cluster_scMethrix("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    
    expect_error(cluster_scMethrix(scm="not scMethrix"),msg.validateExp)
    expect_error(cluster_scMethrix(scm,assay="not an assay"),msg.validateAssay)
    
    scm <- transform_assay(scm,new_assay="fill", trans = function(x) fill(x,fill=0))

    expect_error(cluster_scMethrix(scm,assay="fill",type="not a type"), msg.validateArg)
    
    expect_error(cluster_scMethrix(scm,assay="fill",dist="not a dist"), msg.validateType)
    
    x <- matrix(rnorm(100), nrow = ncol(scm))
    expect_error(cluster_scMethrix(scm,assay="fill",dist=dist(x)), "Invalid distance matrix")
    
    if (is_h5(scm)) {
      expect_warning(dist <- get_distance_matrix(scm, assay="fill"),"Distance matrix cannot be generated for HDF5 data")
    } else {
      dist <- get_distance_matrix(scm, assay="fill")
    }
    
    invisible(lapply(types, function(type) {
      scm.c <- scm
      if (type == "model") {
        expect_warning(
          scm.c <- cluster_scMethrix(scm.c, dist = dist, n_clusters=n_clusters,assay="fill",type=type,colname=name)
          ,"n_clusters is ignored")
      } else {
        scm.c <- cluster_scMethrix(scm.c, dist = dist, n_clusters=n_clusters,assay="fill",type=type,colname=name)
      }
      
      cd <- colData(scm.c)
      expect_is(scm.c,"scMethrix")
      expect_equivalent(colnames(cd),c(colnames(colData(scm)),name))
      expect_true(all(cd$Cluster %in% 1:n_clusters))
    }))
  }))
})

test_that("append_colData", {

  expect_error(append_colData("not scMethrix"),msg.validateExp)
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
