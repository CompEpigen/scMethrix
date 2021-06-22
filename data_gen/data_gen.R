
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
files <- list.files (dir,full.names = TRUE)
files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]

m <- read_beds(files=files,h5=FALSE,cov_idx=c(5,6))
m <- subset_scMethrix(m,contigs=c("chr1","chr2"))
m <- mask_scMethrix(m, low_count = NULL, high_quantile=0.9)

cpg <- lapply(list(0,1,2,3), function(n) {

  set.seed(123)
  c <- which(rowSums(is.na(get_matrix(m)))==n)
  c <- sample(c,50)
  
})
 
m <- m[sort(unlist(cpg))]
suffix <- "_sub"

export_bed(m, path = dir, suffix=suffix)

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
files <- list.files (dir,full.names = TRUE)
files <- files[grepl(paste0(".*",suffix,".bedgraph$"), files,ignore.case = TRUE)]

m <- read_beds(files=files,h5=FALSE,cov_idx=5)