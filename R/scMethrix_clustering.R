#--- get_distance_matrix -------------------------------------------------------------------------------------
#' Get the pair-wise distance matrix for an assay
#' @details Utilizes mainly the bioDist package to determine various distance metrics to be used for later clustering. 
#' @param scm scMethrix; Input \code{\link{scMethrix}} object
#' @param assay string; The assay to use. Default is 'score'
#' @param type string; The type of distance metric to use. Available options are "pearson", "spearman", "tau", "euclidean", "manhattan", "canberra", "binary", "minkowski". An aribitrary distance function can also be used, so long as the input takes just the specified matrix.
#' @param verbose boolean; flag to output messages or not
#' @return matrix; the distance matrix
#' @seealso <https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/dist>
#' @seealso <https://www.bioconductor.org/packages//2.7/bioc/html/bioDist.html>
#' @examples
#' data('scMethrix_data')
#' scMethrix_data <- impute_regions(scMethrix_data) 
#' get_distance_matrix(scMethrix_data,assay = "impute") 
#' 
#' # For an arbitrary distance function
#' library(bioDist)
#' fun <- bioDist::spearman.dist
#' get_distance_matrix(scMethrix_data,assay = "impute",type = bioDist::spearman.dist)
#' @export
get_distance_matrix <- function(scm, assay="score",type=c("pearson", "spearman", "kendall", "euclidean", "manhattan", "canberra", "binary", "minkowski"),verbose=TRUE) {
  
  if (!requireNamespace("bioDist", quietly = TRUE)) {
    stop("Package \"bioDist\" must be installed to use this function.", call. = FALSE)
  }
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  if (!is.function(type)) type = .validateArg(type,get_distance_matrix)
  .validateType(verbose,"boolean")
  
  if (any(is.na(get_matrix(scm,assay=assay)))) stop("There are NA values present in the matrix. Please fill/impute/bin to remove NAs.")
  
  if (is_h5(scm)) {
    warning("Distance matrix cannot be generated for HDF5 data. Data will be cast as matrix for imputation.")
  }
  
  #---- Function code ------------------------------------------------------
  mtx <- as.matrix(t(get_matrix(scm,assay=assay)))
  
  if (is.function(type)) { # For the arbitrary case
    dist <- type(mtx)
  } else if (type == "spearman") {
    dist <- bioDist::spearman.dist(mtx)
  } else if (type == "pearson") {
    dist <- bioDist::cor.dist(mtx)
  } else if (type == "kendall") {
    dist <- bioDist::tau.dist(mtx)
  } else if (type %in% c("euclidean", "manhattan", "canberra", "binary", "minkowski")) {
    dist <- stats::dist(mtx, method = type)
  } else {
    stop("Invalid distance metric specified")
  }
  
  if (all(dim(as.matrix(dist)) != rep(ncol(scm),2))) {
    stop("Outputted distance matrix size does not match number of samples. This should not happen.")
  }
  
  return(dist)
}

#---- cluster_scMethrix ------------------------------------------------------------------------------------------------
#' Generates a cluster object for an \code{\link{scMethrix}} object
#' @details Enables multiple methods of clustering to classify samples in an \code{\link{scMethrix}} object. Either an \code{\link{scMethrix}} object or a \code{\link[stats]{dist}} object must be provided for clustering.
#' @param scm scMethrix; Input \code{\link{scMethrix}} object. If this is specified the distance matrix will be a generic \code{\link[bioDist]{spearman.dist}} distance
#' @param dist dist; Optional. A distance matrix generated for an assay. Will use default paramaters for \code{\link{get_distance_matrix}}.
#' @param assay string; The assay to use. Default is 'score'
#' @param type string; The type of distance metric to use. Available options are 'hierarchical', 'partition', "model". An arbitrary cluster function can be used, and must return a named vector containing integers representing the cluster membership (e.g. \code{c(C1=1,C2=1,C3=1,C4=2)}).
#' @param n_clusters integer; the desired number of clusters. This is ignored for model-based clustering
#' @param colname string; the name of the colData column that contains the cluster information
#' @param verbose boolean; flag to output messages or not
#' @param ... Additional parameters for the clustering functions
#' @return An \code{\link{scMethrix}} object
#' @seealso [get_distance_matrix()] for distance metrics, [hclust()] for heirarchical clustering, [kmeans()] for partition clustering, [mclust::Mclust()] for model clustering 
#' @examples
#' data('scMethrix_data')
#' scMethrix_data <- impute_regions(scMethrix_data)
#' dist <- get_distance_matrix(scMethrix_data,assay = "impute")  
#' 
#' # For a generic clustering function
#' # The function must return a named vector of integers
#' fun <- function (dist) {
#'     fit <- hclust(dist, method="ward.D")
#'     fit <- cutree(fit, k=2)
#'     return(fit)
#'     }
#' 
#' fun(dist) # Example of arbitrary function output
#' colData(cluster_scMethrix(scMethrix_data, dist = dist, type = fun))
#' @export
cluster_scMethrix <- function(scm = NULL, dist = NULL,  assay="score", type=c("hierarchical", "partition", "model"),
                              colname = "Cluster", n_clusters = 1, verbose = TRUE, ...) {

  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  if (typeof(type) != "closure") {
    .validateType(type,"String")
    type = .validateArg(type,cluster_scMethrix)
  }
  .validateType(n_clusters,"integer")
  .validateType(colname,"string")
  .validateType(verbose,"boolean")
  .validateType(dist,c("dist","null"))

  if (is.null(dist)) {
    dist <- get_distance_matrix(scm, assay=assay)
  } else {
    if (attr(dist,"Size") != ncol(scm) || !setequal(labels(dist),sampleNames(scm))) 
      stop("Invalid distance matrix. Must contain all samples present in the experiment")
  }
  
  Cluster <- Sample <- NULL
  
  #---- Function code ------------------------------------------------------

  if (.validateType(type,"function",throws=F)) {
    
    fit <- type(dist)
    if (!setequal(labels(fit),sampleNames(scm))) stop("Invalid cluster function. Must output a named vector containing all samples in the experiment.")
    colData <- data.frame(Sample = names(fit), Cluster = fit)
  } else if (type=="hierarchical") {
    fit <- stats::hclust(dist, ...)
    fit <- stats::cutree(fit, k=n_clusters)
    colData <- data.frame(Sample = names(fit), Cluster = fit)
  } else if (type=="partition") {
    fit <- stats::kmeans(dist, centers = min(n_clusters,attr(dist,"Size")-1)) # Max clusters = nrow(scm)-1
    colData <- data.frame(Sample = names(fit$cluster), Cluster = fit$cluster)
  } else if (type == "model") {
    
    if (!requireNamespace("mclust", quietly = TRUE)) {
      stop("Package \"mclust\" must be installed to use this function.", call. = FALSE)
    }
    
    if (!is.null(n_clusters)) warning("n_clusters is ignored for model-based clustering")
    fit <- mclust::Mclust(as.matrix(t(get_matrix(scm,assay=assay))), ...)
    colData <- data.frame(Sample = names(fit$classification), Cluster = fit$classification)
  }
  # else if (type == "density") {
  #   if (!is.null(n_clusters)) warning("n_clusters is ignored for density clustering")
  #   fit <- dbscan(dist, eps, minPts = n_clusters, borderPoints = TRUE, ...)
  #   
  # }
  # 
  
  row.names(colData) <- colData$Sample
  colData <- subset(colData, select=-Sample)

  names(colData)[names(colData) == "Cluster"] <- colname
  
  scm <- append_colData(scm,colData)
  
  validObject(scm)
  return(scm)
}

#---- append_colData ---------------------------------------------------------------------------------------------------
#' Appends colData in an scMethrix object
#' @details Typically used for clustering. Allows additional information to be added to colData in an scMethrix object after the object creation. It does this via a left join on the original colData. Any samples not included in the colData object will be filled with NAs.
#' @param scm scMethrix; Input \code{\link{scMethrix}} object
#' @param colData dataframe-like or named vector; For a dataframe-like, must contain row names that correspond with the input \code{\link{scMethrix}} object. For a named vector, vector names must correspond to the row names of the input \code{\link{scMethrix}} object
#' @param name string; the name of the column for named vector input. Ignored for matrix input
#' @return An \code{\link{scMethrix}} object
#' @examples
#' 
#' # For dataframe-like input
#' data('scMethrix_data')
#' colData <- colData(scMethrix_data)
#' colData["Cluster"] <- 1:nrow(colData)
#' scMethrix_data <- append_colData(scMethrix_data,colData)
#' colData(scMethrix_data)
#' 
#' # For named vector input
#' data('scMethrix_data')
#' colData <- c(C1=1,C2=1,C3=1,C4=2)
#' scMethrix_data <- append_colData(scMethrix_data,colData, name="Cluster")
#' colData(scMethrix_data)
#' 
#' @export
append_colData <- function(scm = NULL, colData = NULL, name = "Data") {

  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  #.validateType(colData,c("vector","dataframe","S4")) #TODO: fix this or remove functino totally
  .validateType(name,"string")
  
  Row.names <- NULL
  #---- Function code ------------------------------------------------------
  # Convert vector to data.frame
  if (is.vector(colData)) {
    colData <- as.data.frame(colData)
    colnames(colData) <- name
  }

  cd <- colData(scm)
  
  cols <- intersect(colnames(colData), colnames(colData(scm)))

  if (length(cols) > 0) {
    warning("Colnames of colData already exist in scMethrix object (",paste(cols,collapse=", "),"). These will be overwritten.")
    cd <- cd[, !(colnames(cd) %in% cols), drop=FALSE]
  }
  
  n_samples <- length(intersect(row.names(cd),row.names(colData)))
  
  if (n_samples != nrow(cd)) warning(nrow(cd)-n_samples," samples are not specified in colData")

  cd <- merge(cd,colData,by=0, all.x = TRUE)
  row.names(cd) <- cd$Row.names
  cd <- subset(cd, select=-c(Row.names))
  colData(scm) <- cd
  
  validObject(scm)
  return(scm)
}
