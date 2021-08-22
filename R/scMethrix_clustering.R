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
#' @seealso https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/dist
#' @seealso https://www.bioconductor.org/packages//2.7/bioc/html/bioDist.html
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

