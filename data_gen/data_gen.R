
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
files <- list.files (dir,full.names = TRUE)
files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]

m <- read_beds(files=files,h5=FALSE,cov_idx=c(5,6))

cpg <- lapply(list(0,1,2,3), function(n) {

  set.seed(123)
  c <- which(rowSums(is.na(get_matrix(m)))==n)
  c <- sample(c,25)
  
})
 
m <- m[sort(unlist(cpg))]
suffix <- "_sub"

export_bed(m = m, path = dir,suffix=suffix)

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
files <- list.files (dir,full.names = TRUE)
files <- files[grepl(paste0(".*",suffix,".bedgraph$"), files,ignore.case = TRUE)]

m <- read_beds(files=files,h5=FALSE,coverage=TRUE)