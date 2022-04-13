# Data is randomly generated from Luo, C. et al. (2017). Single-cell methylomes identify neuronal subtypes and regulatory 
# Using files: final_SRR5388959.sra_1_CpG.bedgraph, final_SRR5389083.sra_1_CpG.bedgraph,
# final_SRR5389125.sra_1_CpG.bedgraph, final_SRR5389145.sra_1_CpG.bedgraph

set.seed(123)
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
src_files <- list.files (dir,full.names = TRUE)
src_files <- src_files[grepl(".*bedgraph$", src_files,ignore.case = TRUE)]

library(BSgenome.Mmusculus.UCSC.mm10)
ref_cpgs <- suppressWarnings(scMethrix::extract_CpGs(genome = "BSgenome.Mmusculus.UCSC.mm10"))

#read <- read_beds(files=src_files,h5=FALSE,cov_idx=c(5,6))
reads <- read_beds(files=src_files, chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, M_idx = 5, U_idx=6, ref_cpgs = ref_cpgs)
reads <- subset_scMethrix(reads,contigs=paste0("chr",1:5))
#m <- mask_scMethrix(m, low_count = NULL, high_quantile=0.9)

cpg <- lapply(0:length(src_files), function(n) {
  c <- which(rowSums(is.na(counts(reads)))==n)
  c <- sample(c,sample(50:100, 1))
})
 
reads <- reads[sort(unlist(cpg))]
suffix <- "_sub"

exportBed(reads, path = dir, assays = c("score","counts"), suffix=suffix, na.rm = FALSE)

files <- list.files (dir,full.names = TRUE)
files <- files[grepl(paste0(".*.bed$"), files,ignore.case = TRUE)]

file.rename(files, paste0(dir,"/C", 1:length(files),".bed"))

files <- list.files (dir,full.names = TRUE)
files <- files[grepl(".*bed$", files,ignore.case = TRUE)]
files <- setdiff(files, src_files)

colData <- data.frame(phenotype = c(rep("C",2),rep("N",2)))
row.names(colData) <- paste0("C",1:length(files))

scm <- read_beds(files=files, chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)

#m <- read_beds(files=files,h5=FALSE,cov_idx=5,colData = colData)