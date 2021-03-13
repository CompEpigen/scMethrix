list.of.packages <- c("data.table","plyr","HDF5Array","tictoc","beepr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

setwd("D:\\Documents\\School\\Thesis\\scMethrix\\sample.data\\bedtools")

setwd("D:/Documents/School/Thesis/scMethrix/sample.data/Yunhee/GSE97179")


files <- list.files (getwd(),full.names = TRUE)
files <- files[1:100]