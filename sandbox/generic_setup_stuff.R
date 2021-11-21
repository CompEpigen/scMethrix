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

#Psuedo Analysis
library(AnnotationHub)
ah = AnnotationHub()
#qhs = query(ah, c("RefSeq", "Mus musculus", "mm10"))
qhs = query(ah, c("RefSeq", "Homo sapiens", "hg19"))
genes = qhs[[1]]
proms = promoters(genes, upstream=1000, downstream=1000)

hg38 <- ah[["AH97949"]]
proms = promoters(genes(hg38), upstream=1000, downstream=1000)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes(hg38)
proms <- promoters(hg38)

reg <- reduce(c(genes(hg38),promoters(hg38)))

scm <- cbind(scm1,scm2)
scm <- cbind(scm,scm3)

scm.genic <- subset_scMethrix(scm,region = reg)
saveHDF5SummarizedExperiment(scm.genic,dir="D:/Git/sampleData/Gaiti/bin.genic")


scm.prom <- bin_scMethrix(scm,regions = promoters(hg38),bin_size=10000000,h5_dir = "D:/Git/sampleData/Gaiti/bin.prom")
scm.gene <- bin_scMethrix(scm,regions = genes(hg38),bin_size=10000000,h5_dir = "D:/Git/sampleData/Gaiti/bin.gene")

scm <- load_scMethrix(dir="D:/Git/sampleData/3samp.h5")
scm <- mask_by_coverage(scm,avg_threshold=2)
scm <- bin_scMethrix(scm,proms,h5_dir = paste0(tempdir(),"/h5bin"))
scm <- convert_HDF5_scMethrix(scm)
saveRDS(scm, file = "D:/Git/sampleData/workingDir/scm_bin_100k.rds")
scm.bin <- readRDS("D:/Git/sampleData/workingDir/scm_bin_100k.rds")
scm.bin <- mask_by_sample(scm.bin,prop_threshold=.95)
scm.bin <- remove_uncovered(scm.bin)
scm.impute <- impute_regions(scm.bin,type="kNN")
scm.binarize <- transform_assay(scm.impute,trans=binarize,assay="impute",new_assay="bin")
scm.umap <- dim_red_scMethrix(scm.impute, assay="impute",type="tSNE",top_var = nrow(scm.impute))

CairoWin(width=5,height=5)
plot_dim_red(scm.umap,"tSNE",col_anno = "Cell") + ggtitle("tSNE") + geom_point(alpha = 1/10) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=20))
dist <- get_distance_matrix(scm.umap,assay="impute",type="euclidean")

scm.cluster <- scm.umap

scm.cluster <- cluster_scMethrix(scm.umap, dist = dist, 
               assay = "impute",type = "model",colname = "model")
plot_dim_red(scm.cluster,"tSNE",col_anno = "model") + ggtitle("Model-based") + geom_point(alpha = 1/10) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=20)) +
  scale_fill_discrete(name = "Cluster")

scm.cluster <- cluster_scMethrix(scm.cluster,dist = dist, assay="impute",type="heir",colname="heir",n_clusters=4)
plot_dim_red(scm.cluster,"tSNE",col_anno = "heir") + ggtitle("Heirarchical") + geom_point(alpha = 1/10) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=20))

scm.cluster <- cluster_scMethrix(scm.cluster,dist = dist, assay="impute",type="part",colname="part",n_clusters=3)
plot_dim_red(scm.cluster,"tSNE",col_anno = "part") + ggtitle("Partitioned") + geom_point(alpha = 1/10) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=20))
