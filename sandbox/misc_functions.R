library(BRGenomics)
library(data.table)
library(GenomicRanges)
library(tictoc)
library(beepr)
library(SingleCellExperiment)
contigs <- c("chr1","chr3")
regions <- GRanges(seqnames = "chr1", ranges = IRanges(1,150)) 
samples <- c("bed1","bed3")

tic();subset_scMethrix(tbx.scm,contigs = contigs);toc();beep()

# Load locations
setwd("D:\\Documents\\School\\Thesis\\scMethrix\\sample.data\\s50cell\\bedgraph")


install.packages("rbenchmark")

library(seqminer)
library(rbenchmark)

benchmark(
  "tabix" = {
    
    tbx <- Rsamtools::scanTabix(file)
   
    tbx <- data.frame(lapply(tbx,function(x) {str_split_fixed(x, "\t", 4)}))
    
    colnames(tbx) <- c("chr","start","end",strsplit(basename(file), "[.]")[[1]][1])
    tbx[,2] <- as.integer(tbx[,2])
    tbx[,3] <- as.integer(tbx[,3])
    tbx[,4] <- as.integer(tbx[,4])
    
  },
  "seqminer" = {
   # seq <- seqminer::tabix.read.table(file, col.names = FALSE)
    
  },
  
  "fread" = {
    fre <-
      data.table::fread(file, header = FALSE)
    
    fre <- data.frame(fre)
    
    colnames(fre) <- c("chr","start","end",strsplit(basename(file), "[.]")[[1]][1])
  
  },
  replications = 1,
  columns = c(
    "test",
    "replications",
    "elapsed",
    "relative",
    "user.self",
    "sys.self"
  )
)

beep()