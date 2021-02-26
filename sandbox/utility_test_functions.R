### Generic functions to  help with testing #################################
sp2mat <- function(m) {
  
  return(as.matrix(m))  
  
  
}

mat2sp <- function(m) {
  
  return(Matrix(m, sparse=TRUE))
  
}
