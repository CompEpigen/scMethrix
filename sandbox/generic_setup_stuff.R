list.of.packages <- c("SingleCellExperiment","data.table","plyr","HDF5Array","tictoc","beepr",
                      "GenomicRanges","parallel","roxygen2","dplyr","rbenchmark","testthat","rtracklayer",
                      "tools","microbenchmark","measurements","magrittr","doParallel","parallel",
                      "Cairo","ggplot2","methrix","BSgenome","BSgenome.Hsapiens.UCSC.hg19","usethis",
                      "BSgenome.Mmusculus.UCSC.mm10","pkgdown","umap","stringi","missMDA","Rtsne","missForest",
                      "impute")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages,new.packages)

source("D:/Git/scMethrix/R/accessory_funcs.R")
source("D:/Git/scMethrix/R/scMethrix_operations.R")
source("D:/Git/scMethrix/R/scMethrix_object.R")
source("D:/Git/scMethrix/R/read_beds.R")
source("D:/Git/scMethrix/R/scMethrix_dimensionality.R")
source("D:/Git/scMethrix/R/scMethrix_transforms.R")
source("D:/Git/scMethrix/R/scMethrix_plot.R")
source("D:/Git/scMethrix/tests/testthat/setup.R")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

#setwd("D:/Documents/School/Thesis/scMethrix/sample.data/small/")

## Reference CpG sets
Hg19_cpgs <- methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
mm10_cpgs <- methrix::extract_CPGs(ref_genome = "BSgenome.Mmusculus.UCSC.mm10")
mm10_cpgs <- mm10_cpgs$cpgs[,1:3]
ref_cpgs <- mm10_cpgs

# Generic data import
setwd("D:/Git/sampleData/Yunhee.GSE97179")

setwd("D:/Git/sampleData/test_data")

setwd("D:/Git/sampleData/100cell")

files <- list.files (getwd(),full.names = TRUE)

files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]

files <- files[1:20]

scm.big.h5 <- read_beds(files=files,h5=TRUE,h5_dir=paste0(getwd(),"/sse"),cov=c(5,6),replace=TRUE)
scm.big.mem <- read_beds(files=files,h5=FALSE)

scm.big.mem <- read_beds(files=files,h5=FALSE, ref_cpgs = mm10_cpgs$cpgs, n_threads = 8)

scm.big.h5 <- read_beds(files=files,h5=TRUE,h5_dir=paste0(getwd(),"/sse"),cov=c(5),
                        replace=TRUE, ref_cpgs = mm10_cpgs$cpgs, n_threads = 8)

scm.20.mem <- readRDS(file = "D:/Git/sampleData/scm.20.mem.rds")

# Useful functions
devtools::test()
devtools::install()
pkgdown::build_site()