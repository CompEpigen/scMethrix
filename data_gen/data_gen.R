
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
files <- list.files (dir,full.names = TRUE)
files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]

m <- read_beds(files=files,h5=FALSE,coverage=TRUE)

cpg <- lapply(c(0,1,2,3), function(n) {
  
  set.seed(123)
  c <- which(rowSums(is.na(get_matrix(m)))==n)
  c <- sample(c,25)
  
})

cpg <- sort(unlist(cpg))
data <- m[cpg]

export_bed(m = data,path = dir,suffix="_sub")