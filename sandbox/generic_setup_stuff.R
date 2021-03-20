list.of.packages <- c("SingleCellExperiment","data.table","plyr","HDF5Array","tictoc","beepr",
                      "GenomicRanges","parallel","roxygen2","dplyr","rbenchmark","testthat","rtracklayer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages,new.packages)

setwd("D:\\Documents\\School\\Thesis\\scMethrix\\sample.data\\bedtools")

setwd("D:/Documents/School/Thesis/scMethrix/sample.data/Yunhee/GSE97179")

setwd("D:/Documents/School/Thesis/scMethrix/sample.data/small/")

files <- list.files (getwd(),full.names = TRUE)

files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]

files <- files[1:10]

h5_dir <- getwd()

tic()
  
scm <- read_beds(files=files,h5=TRUE,h5_dir=paste0(getwd(),"/sse"))

scm <- read_beds(files=files,h5=FALSE)
  
toc()

