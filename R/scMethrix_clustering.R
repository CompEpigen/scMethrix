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
# groups <- cutree(fit, k=5) # cut tree into 5 clusters
# # draw dendogram with red borders around the 5 clusters
# rect.hclust(fit, k=5, border="red")
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