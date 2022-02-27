list.of.packages <- c("SingleCellExperiment","data.table","plyr","HDF5Array","tictoc","beepr",
                      "GenomicRanges","parallel","roxygen2","dplyr","rbenchmark","testthat","rtracklayer",
                      "tools","microbenchmark","measurements","magrittr","doParallel","parallel",
                      "Cairo","ggplot2","methrix","BSgenome","usethis",
                      "pkgdown","umap","stringi","missMDA","Rtsne","missForest",
                      "impute","profvis",'Melissa','Metrics','SimDesign','bioDist','dbscan','AnnotationHub','mclust',
                      'bsseq','Seurat')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
  BiocManager::install(new.packages)
}
status <- lapply(list.of.packages, require, character.only = TRUE)
names(status) <- list.of.packages
suppressWarnings(if (!all(status)) status[which(status==FALSE)])
rm(list.of.packages,new.packages)

# assign("time.all", numeric(), envir=topenv())
# assign("time.split", numeric(), envir=topenv())

if (!exists("ref_cpgs")) ref_cpgs2 <- readRDS("D:/Git/sampleData/ref_cpgs.rds")

source("D:/Git/scMethrix/R/zzz.R")
source("D:/Git/scMethrix/R/accessory_funcs.R")
source("D:/Git/scMethrix/R/scMethrix_operations.R")
source("D:/Git/scMethrix/R/scMethrix_object.R")
source("D:/Git/scMethrix/R/read_beds.R")
source("D:/Git/scMethrix/R/scMethrix_dimensionality.R")
source("D:/Git/scMethrix/R/scMethrix_transforms.R")
source("D:/Git/scMethrix/R/scMethrix_clustering.R")
source("D:/Git/scMethrix/R/scMethrix_plot.R")
source("D:/Git/scMethrix/R/scMethrix_export.R")
source("D:/Git/scMethrix/tests/testthat/setup.R")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

#setwd("D:/Documents/School/Thesis/scMethrix/sample.data/small/")

## Reference CpG sets
Hg19_cpgs <- extract_CpGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
Hg38_cpgs <- extract_CpGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
mm10_cpgs <- extract_CpGs(ref_genome = "BSgenome.Mmusculus.UCSC.mm10")
ref_cpgs <- mm10_cpgs
rm(mm10_cpgs)
saveRDS(ref_cpgs, file = "D:/Git/sampleData/ref_cpgs.rds")

# Generic data import
setwd("D:/Git/sampleData/Yunhee.GSE97179")

setwd("D:/Git/sampleData/test_data")

setwd("D:/Git/sampleData/100cell")

setwd("D:/Git/sampleData/mini")

setwd("D:/Git/sampleData/3samp")

files <- list.files (getwd(),full.names = TRUE)

files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]

files <- files[1:590]

col_list <- parse_source_idx(chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, M_idx=5, U_idx=6)

#With Coverage
scm.big.h5 <- read_beds(files=files,h5=TRUE,h5_dir=paste0(tempdir(),"/sse"),ref_cpgs = ref_cpgs, replace=TRUE,
                        chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, M_idx=5, U_idx=6, colData = colData, n_threads=0, batch_size = 10)
saveHDF5SummarizedExperiment(scm.big.h5,dir="D:/Git/sampleData/3samp.h5",replace=TRUE)
scm.big.h5 <- load_scMethrix(dir="D:/Git/sampleData/500.h5")
scm.big.h5 <- load_scMethrix(dir="D:/Git/sampleData/3samp.h5")

scm.big.mem <- read_beds(files=files,h5=FALSE,n_threads = 1, colData = colData,
                         chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, M_idx=5, U_idx=6)

#Without coverage

scm.big.mem <- read_beds(files=files,h5=FALSE, ref_cpgs = mm10_cpgs$cpgs, n_threads = 8, 
                         chr_idx=1, start_idx=2, end_idx=3, beta_idx=4)
scm.big.h5 <- read_beds(files=files,h5=TRUE, h5_dir=paste0(tempdir(),"/sse"), replace=TRUE, ref_cpgs = ref_cpgs,
                        chr_idx=1, start_idx=2, end_idx=3, beta_idx=4,
                        colData = colData, n_threads=0, batch_size = 10)

#Methrix input
setwd("D:/Documents/School/Thesis/methrix/methrix_data_generation")
scm.methrix <- read_beds(files=files,h5=FALSE, chr_idx=1, start_idx=2, M_idx=3, U_idx=4)

# Useful functions
devtools::test()
devtools::install()
devtools::run_examples()
devtools::load_all(".")

setwd("D:/Git/scMethrix")
pkgdown::build_site()

saveRDS(colData, file = "colData.rds")

#Testing stuff
scm1 <- load_scMethrix("d:/scm1.rds")
scm2 <- load_scMethrix("d:/scm2.rds")
rr1 = rowRanges(scm1); rr2 = rowRanges(scm2)
scm <- merge_scMethrix2(scm1,scm2,by_row_name = T)


