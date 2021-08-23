# # Ward Hierarchical Clustering
# tic()
# 
# d <- dist(score(scm.impute), method = "euclidean") # distance matrix
# toc()
# 
# tic()
# fit <- hclust(d, method="ward.D")
# toc()
# plot(fit) # display dendogram
# groups <- cutree(fit, k=4) # cut tree into 5 clusters
# # draw dendogram with red borders around the 5 clusters
# rect.hclust(fit, k=2, border="red")
# 
# 
# 
# # Model Based Clustering
# library(mclust)
# tic()
# fit <- Mclust(get_matrix(scm.dim,assay="impute"))
# toc()
# plot(fit) # plot results
# summary(fit) # display the best model

#' Get the pair-wise distance matrix for an assay
#' @details Utilizes mainly the bioDist package to determine various distance metrics to be used for later clustering. 
#' @param scm scMethrix; Input \code{\link{scMethrix}} object
#' @param assay string; The assay to use. Default is 'score'
#' @param type string; The type of distance metric to use. Available options are 'spearman' and 'tau'. 
#' @param verbose boolean; flag to output messages or not
#' @return matrix; the distance matrix
#' @import bioDist
#' @seealso <https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/dist>
#' @seealso <https://www.bioconductor.org/packages//2.7/bioc/html/bioDist.html>
#' @examples
#' @export
get_distance_matrix <- function(scm, assay="score",type="spearman",verbose=TRUE) {
  
  # mtx <- as.data.table(get_matrix(scm,assay=assay))
  # mtx <- setDT(transpose(mtx,keep.names="sample"))[]
  # mtx <- as.matrix(mtx)
  # 
  
  mtx <- as.matrix(t(get_matrix(scm,assay=assay)))
  
  if (any(is.na(mtx))) stop("There are NA values present in the matrix. Please fill/impute/bin to remove NAs.")
  
  if (type == "spearman") {
    return (spearman.dist(mtx))
  } else if (type == "tau") {
    return(tau.dist(mtx))
  } 
  
}

#' Generates a cluster object for an \code{\link{scMethrix}} object
#' @details Enables multiple methods of clustering to classify samples in an \code{\link{scMethrix}} object. Either an \code{\link{scMethrix}} object or a \code{\link[stats]{dist}} object must be provided for clustering.
#' @param scm scMethrix; Input \code{\link{scMethrix}} object. If this is specified the distance matrix will be a generic \code{\link[bioDist]{spearman.dist}} distance
#' @param dist dist; A distance matrix generated for an scMethrix object. Not necessary
#' @param assay string; The assay to use. Default is 'score'
#' @param type string; The type of distance metric to use. Available options are 'spearman' and 'tau'. 
#' @param verbose boolean; flag to output messages or not
#' @return matrix; the distance matrix
#' @import bioDist
#' @seealso [get_distance_matrix()]
#' @examples
#' @export
cluster_scMethrix <- function(scm = NULL, dist = NULL, assay="score", verbose = TRUE, type="hierarchical") {

  if (is.null(dist)) {
    if (is.null(scm)) {
      stop("Either scm or dist must be specified")
    } else {
      dist <- get_distance_matrix(scm, assay=assay)
    }
  }
  
  if (type=="hierarchical") {
    fit <- hclust(dist, method="ward.D")
    return(fit)
  } 
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
  
  if (!("Sample" %in% colnames(colData))) {
    colData["Sample"] <- rownames(colData)
  }
  
  cd <- colData(scm)
  cd["Sample"] <- rownames(cd)
  
  n_samples <- length(intersect(cd$Sample,colData$Sample))
  
  if (n_samples != nrow(cd)) warning(nrow(cd)-n_samples," samples are not specified in colData")

  cd <- merge(cd,colData,by="Sample", all.x = TRUE)
  
  row.names(cd) <- cd$Sample
  cd <- within(cd, rm(Sample))
  
  colData(scm) <- cd
  
  return(scm)
}
