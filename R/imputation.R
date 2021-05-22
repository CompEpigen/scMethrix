impute_scMethrix <- function (m) {
  
  impute <- transform_assay(m,assay = "score",name = "impute",trans = function(x) {ifelse(m > 50,1,0)})

  
}






