list.of.packages <- c("SingleCellExperiment","data.table","plyr","HDF5Array","tictoc","beepr",
                      "GenomicRanges","parallel","roxygen2","dplyr","rbenchmark","testthat","rtracklayer",
                      "tools","microbenchmark","measurements","magrittr","doParallel","parallel",
                      "Cairo","ggplot2","methrix","BSgenome","BSgenome.Hsapiens.UCSC.hg19","usethis",
                      "BSgenome.Mmusculus.UCSC.mm10","pkgdown","umap","stringi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages,new.packages)

#setwd("D:/Documents/School/Thesis/scMethrix/sample.data/small/")

## Reference CpG sets
Hg19_cpgs <- methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
mm10_cpgs <- methrix::extract_CPGs(ref_genome = "BSgenome.Mmusculus.UCSC.mm10")
mm10_cpgs <- mm10_cpgs$cpgs[,1:3]
ref_cpgs <- mm10_cpgs

# Generic data import
setwd("D:/Git/sampleData/Yunhee.GSE97179")

setwd("D:/Git/sampleData/test_data")

files <- list.files (getwd(),full.names = TRUE)

files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]

files <- files[1:50]

scm.big.h5 <- read_beds(files=files,h5=TRUE,h5_dir=paste0(getwd(),"/sse"),cov=c(5,6),replace=TRUE)

scm.big.mem <- read_beds(files=files,h5=FALSE)

# Useful functions
devtools::test()
devtools::install()
pkgdown::build_site()