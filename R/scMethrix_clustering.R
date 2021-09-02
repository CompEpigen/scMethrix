#' Get the pair-wise distance matrix for an assay
#' @details Utilizes mainly the bioDist package to determine various distance metrics to be used for later clustering. 
#' @param scm scMethrix; Input \code{\link{scMethrix}} object
#' @param assay string; The assay to use. Default is 'score'
#' @param type string; The type of distance metric to use. Available options are "pearson", spearman", "tau", "euclidean". "maximum", "manhattan", "canberra", "binary", "minkowski". An aribitrary distance function can also be used, so long as the input takes just the specified matrix.
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
  
  mtx <- as.matrix(t(get_matrix(scm,assay=assay)))
  
  if (any(is.na(mtx))) stop("There are NA values present in the matrix. Please fill/impute/bin to remove NAs.")
  
  if (typeof(type) == "closure") { # For the arbitrary case
    dist <- type(mtx)
  } else if (type == "spearman") {
    dist <- spearman.dist(mtx)
  } else if (type == "pearson") {
    dist <- cor.dist(mtx)
  } else if (type == "tau") {
    dist <- tau.dist(mtx)
  } else if (type %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
    dist <- dist(mtx, method = type)
  } else {
    stop("Invalid distance metric specified")
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
cluster_scMethrix <- function(scm = NULL, dist = NULL, n_clusters = NULL, assay="score", verbose = TRUE, type="hierarchical", ...) {

  Cluster <- Sample <- NULL
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is.null(dist)) dist <- get_distance_matrix(scm, assay=assay)
  
  if (is.null(n_clusters)) n_clusters = attr(dist,"Size")
  
  if (typeof(type) == "closure") {
    fit <- type(dist)
    colData <- data.frame(Sample = names(fit), Cluster = fit)
  } else if (type=="hierarchical") {
    fit <- stats::hclust(dist, method="ward.D", ...)
    fit <- stats::cutree(fit, k=n_clusters)
    colData <- data.frame(Sample = names(fit), Cluster = fit)
  } else if (type=="partition") {
    fit <- stats::kmeans(dist, centers = min(n_clusters,attr(dist,"Size")-1)) # Max clusters = nrow(scm)-1
    colData <- data.frame(Sample = names(fit$cluster), Cluster = fit$cluster)
  } else if (type == "model") {
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
  scm <- append_col_data(scm,colData)
  return(scm)
}

#' Appends colData in an scMethrix object
#' @details Typically used for clustering. Allows additional information to be added to colData in an scMethrix object after the object creation. It does this via a left join on the original colData. Any samples not included in the colData object will be filled with NAs.
#' @param scm scMethrix; Input \code{\link{scMethrix}} object
#' @param colData matrix; A matrix containing colData. Must contain either a column labelled "Sample" or row names that correspond with the input \code{\link{scMethrix}} object
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' colData <- colData(scMethrix_data)
#' colData["Type"] <- "Cell"
#' scMethrix_data <- append_col_data(scMethrix_data,colData)
#' colData(scMethrix_data)
#' @export
append_col_data <- function(scm, colData) {
  
  Sample <- NULL
  
  if (!("Sample" %in% colnames(colData))) {
    colData["Sample"] <- rownames(colData)
  }
  
  cd <- colData(scm)
  cd["Sample"] <- rownames(cd)
  
  n_samples <- length(intersect(cd$Sample,colData$Sample))
  
  if (n_samples != nrow(cd)) warning(nrow(cd)-n_samples," samples are not specified in colData")

  cd <- merge(cd,colData,by="Sample", all.x = TRUE)
  
  row.names(cd) <- cd$Sample
  cd <- subset(cd, select=-c(Sample))
  
  colData(scm) <- cd
  
  return(scm)
}