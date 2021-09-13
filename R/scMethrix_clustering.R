#' Get the pair-wise distance matrix for an assay
#' @details Utilizes mainly the bioDist package to determine various distance metrics to be used for later clustering. 
#' @param scm scMethrix; Input \code{\link{scMethrix}} object
#' @param assay string; The assay to use. Default is 'score'
#' @param type string; The type of distance metric to use. Available options are "pearson", "spearman", "tau", "euclidean", "manhattan", "canberra", "binary", "minkowski". An aribitrary distance function can also be used, so long as the input takes just the specified matrix.
#' @param verbose boolean; flag to output messages or not
#' @return matrix; the distance matrix
#' @import bioDist
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
get_distance_matrix <- function(scm, assay="score",type="euclidean",verbose=TRUE) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }

  if (typeof(type) != "closure" && !(type %in% c("pearson", "spearman", "kendall", "euclidean", "manhattan", "canberra", "binary", "minkowski"))) stop("Invalid type of distance calculation")
  
  mtx <- as.matrix(t(get_matrix(scm,assay=assay)))
  
  if (any(is.na(mtx))) stop("There are NA values present in the matrix. Please fill/impute/bin to remove NAs.")
  
  if (typeof(type) == "closure") { # For the arbitrary case
    dist <- type(mtx)
  } else if (type == "spearman") {
    dist <- spearman.dist(mtx)
  } else if (type == "pearson") {
    dist <- cor.dist(mtx)
  } else if (type == "kendall") {
    dist <- tau.dist(mtx)
  } else if (type %in% c("euclidean", "manhattan", "canberra", "binary", "minkowski")) {
    dist <- dist(mtx, method = type)
  } else {
    stop("Invalid distance metric specified")
  }
  
  if (all(dim(as.matrix(dist)) != rep(ncol(scm),2))) {
    stop("Outputted distance matrix size does not match number of samples. This should not happen.")
  }
  
  return(dist)
}

#' Generates a cluster object for an \code{\link{scMethrix}} object
#' @details Enables multiple methods of clustering to classify samples in an \code{\link{scMethrix}} object. Either an \code{\link{scMethrix}} object or a \code{\link[stats]{dist}} object must be provided for clustering.
#' @param scm scMethrix; Input \code{\link{scMethrix}} object. If this is specified the distance matrix will be a generic \code{\link[bioDist]{spearman.dist}} distance
#' @param dist dist; Optional. A distance matrix generated for an assay. Will use default paramaters for \code{\link{get_distance_matrix}}.
#' @param assay string; The assay to use. Default is 'score'
#' @param type string; The type of distance metric to use. Available options are 'hierarchical', 'partition', "model". An arbitrary cluster function can be used, and must return a named vector containing integers representing the cluster membership (e.g. \code{c(C1=1,C2=1,C3=1,C4=2)})
#' @param n_clusters integer; the desired number of clusters. This is ignored for model-based clustering
#' @param colname string; the name of the colData column that contains the cluster information
#' @param verbose boolean; flag to output messages or not
#' @param ... Additional parameters for the clustering functions
#' @return An \code{\link{scMethrix}} object
#' @import bioDist
#' @import mclust
#' @seealso [get_distance_matrix()] for distance metrics, [hclust()] for heirarchical clustering, [kmeans()] for partition clustering, [Mclust()] for model clustering 
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
#' cluster_scMethrix(scMethrix_data, dist = dist, type = fun)
#' @export
cluster_scMethrix <- function(scm = NULL, dist = NULL, n_clusters = NULL, assay="score", colname = "Cluster", verbose = TRUE, type="heir", ...) {

  Cluster <- Sample <- NULL
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }
  
  if (typeof(type) != "closure" && !(type %in% c("heir", "part", "model"))) stop("Invalid type of clustering")
  
  if (is.null(dist)) dist <- get_distance_matrix(scm, assay=assay)
  
  if (is.null(n_clusters)) n_clusters = attr(dist,"Size")

  if (typeof(type) == "closure") {
    fit <- type(dist)
    colData <- data.frame(Sample = names(fit), Cluster = fit)
  } else if (type=="heir") {
    fit <- stats::hclust(dist, method="ward.D", ...)
    fit <- stats::cutree(fit, k=n_clusters)
    colData <- data.frame(Sample = names(fit), Cluster = fit)
  } else if (type=="part") {
    fit <- stats::kmeans(dist, centers = min(n_clusters,attr(dist,"Size")-1)) # Max clusters = nrow(scm)-1
    colData <- data.frame(Sample = names(fit$cluster), Cluster = fit$cluster)
  } else if (type == "model") {
    if (!is.null(n_clusters)) warning("n_clusters is ignored for model-based clustering")
    fit <- mclust::Mclust(as.matrix(t(get_matrix(scm,assay=assay))), ...)
    colData <- data.frame(Sample = names(fit$classification), Cluster = fit$classification)
  }
  
  names(colData)[names(colData) == "Cluster"] <- colname
  
  # else if (type == "density") {
  #   if (!is.null(n_clusters)) warning("n_clusters is ignored for density clustering")
  #   fit <- dbscan(dist, eps, minPts = n_clusters, borderPoints = TRUE, ...)
  #   
  # }
  # 
  scm <- append_colData(scm,colData)
  return(scm)
}

#' Appends colData in an scMethrix object
#' @details Typically used for clustering. Allows additional information to be added to colData in an scMethrix object after the object creation. It does this via a left join on the original colData. Any samples not included in the colData object will be filled with NAs.
#' @param scm scMethrix; Input \code{\link{scMethrix}} object
#' @param colData dataframe-like or named vector; For a dataframe-like, Must contain either a column labelled "Sample" or row names that correspond with the input \code{\link{scMethrix}} object. For a named vector, vector names must correspond to the row names of the input \code{\link{scMethrix}} object
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

  Sample <- NULL
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  colData <- data.frame(colData) #TODO: test better for inputs
  if (ncol(colData) == 1) colnames(colData) = name
  
  # if (!is(colData, "vector")) {
  #   stop("A valid colData object must be supplied (named vector or dataframe-like).", call. = FALSE)
  # }

  cd <- colData(scm)
  cols <- colnames(colData) %in% colnames(cd)
  
  if (any(cols)) {
    warning("Colnames of colData already exist in scMethrix object (",colnames(colData)[cols],"). These will be overwritten.")
    cd <- cd[ , !(names(cd) %in% colnames(colData)[cols])]
  }
  
  if (!("Sample" %in% colnames(colData))) colData["Sample"] <- rownames(colData)
  
  cd["Sample"] <- rownames(cd)
  n_samples <- length(intersect(cd$Sample,colData$Sample))
  
  if (n_samples != nrow(cd)) warning(nrow(cd)-n_samples," samples are not specified in colData")

  cd <- merge(cd,colData,by="Sample", all.x = TRUE)
  row.names(cd) <- cd$Sample
  cd <- subset(cd, select=-c(Sample))
  colData(scm) <- cd
  
  return(scm)
}