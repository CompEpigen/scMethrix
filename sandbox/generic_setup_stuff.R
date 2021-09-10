list.of.packages <- c("SingleCellExperiment","data.table","plyr","HDF5Array","tictoc","beepr",
                      "GenomicRanges","parallel","roxygen2","dplyr","rbenchmark","testthat","rtracklayer",
                      "tools","microbenchmark","measurements","magrittr","doParallel","parallel",
                      "Cairo","ggplot2","methrix","BSgenome","BSgenome.Hsapiens.UCSC.hg19","usethis",
                      "BSgenome.Mmusculus.UCSC.mm10","pkgdown","umap","stringi","missMDA","Rtsne","missForest",
                      "impute","profvis",'Melissa','Metrics','SimDesign','bioDist','dbscan','AnnotationHub')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
  BiocManager::install(new.packages)
}
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages,new.packages)

assign("time.all", numeric(), envir=topenv())
assign("time.split", numeric(), envir=topenv())

if (is.null(ref_cpgs)) ref_cpgs <- readRDS("D:/Git/sampleData/ref_cpgs.rds")

source("D:/Git/scMethrix/R/accessory_funcs.R")
source("D:/Git/scMethrix/R/scMethrix_operations.R")
source("D:/Git/scMethrix/R/scMethrix_object.R")
source("D:/Git/scMethrix/R/read_beds.R")
source("D:/Git/scMethrix/R/scMethrix_dimensionality.R")
source("D:/Git/scMethrix/R/scMethrix_transforms.R")
source("D:/Git/scMethrix/R/scMethrix_clustering.R")
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
rm(mm10_cpgs)

# Generic data import
setwd("D:/Git/sampleData/Yunhee.GSE97179")

setwd("D:/Git/sampleData/test_data")

setwd("D:/Git/sampleData/100cell")

setwd("D:/Git/sampleData/mini")

files <- list.files (getwd(),full.names = TRUE)

files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]

files <- files[1:20]

col_list <- parse_source_idx(chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, M_idx=5, U_idx=6)

#With Coverage
scm.big.h5 <- read_beds(files=files,h5=TRUE,h5_dir=paste0(tempdir(),"/sse"),ref_cpgs = ref_cpgs, replace=TRUE,
                        chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, M_idx=5, U_idx=6, colData = colData, n_threads=0, batch_size = 20)
scm.big.h5 <- load_HDF5_scMethrix(dir="D:/Git/sampleData/500.h5")
saveHDF5SummarizedExperiment(scm.big.h5,dir="D:/Git/sampleData/500.h5",replace=TRUE)

scm.big.mem <- read_beds(files=files,h5=FALSE,n_threads = 0, colData = colData,
                         chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, M_idx=5, U_idx=6)

#Without coverage

scm.big.mem <- read_beds(files=files,h5=FALSE, ref_cpgs = mm10_cpgs$cpgs, n_threads = 8, 
                         chr_idx=1, start_idx=2, end_idx=3, beta_idx=4)
scm.big.h5 <- read_beds(files=files,h5=TRUE,h5_dir=paste0(getwd(),"/sse"),cov=c(5),
                        replace=TRUE, ref_cpgs = mm10_cpgs$cpgs, n_threads = 8,
                        chr_idx=1, start_idx=2, end_idx=3, beta_idx=4)

#Methrix input
setwd("D:/Documents/School/Thesis/methrix/methrix_data_generation")
scm.methrix <- read_beds(files=files,h5=FALSE, chr_idx=1, start_idx=2, M_idx=3, U_idx=4)

# Useful functions
devtools::test()
devtools::install()
devtools::run_examples()

setwd("D:/Git/scMethrix")
pkgdown::build_site()

saveRDS(colData, file = "colData.rds")

#Psuedo Analysis
ah = AnnotationHub()
qhs = query(ah, c("RefSeq", "Mus musculus", "mm10"))
genes = qhs[[1]]
proms = promoters(genes)

scm <- load_HDF5_scMethrix(dir="D:/Git/sampleData/500.h5")
scm <- scm[,1:25]
scm <- bin_scMethrix(scm,proms,h5_dir = paste0(tempdir(),"/h5"))
scm <- convert_HDF5_scMethrix(scm)
scm <- mask_by_sample(scm,low_threshold=NULL,prop_threshold=0.95)
scm <- mask_by_variance(scm,low_threshold=0.05)
scm <- remove_uncovered(scm)
scm <- impute_regions(scm)
scm <- transform_assay(scm,trans=binarize,assay="impute",new_assay="bin")
scm <- dim_red_scMethrix(scm,assay="impute",type="UMAP")
plot_dim_red(scm,"UMAP")
