# Data is randomly generated from Luo, C. et al. (2017). Single-cell methylomes identify neuronal subtypes and regulatory 
# Using files: final_SRR5388959.sra_1_CpG.bedgraph, final_SRR5389083.sra_1_CpG.bedgraph,
# final_SRR5389125.sra_1_CpG.bedgraph, final_SRR5389145.sra_1_CpG.bedgraph

set.seed(123)
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
src_files <- list.files (dir,full.names = TRUE)
src_files <- src_files[grepl(".*bedgraph$", src_files,ignore.case = TRUE)]

read <- read_beds(files=src_files,h5=FALSE,cov_idx=c(5,6))
m <- subset_scMethrix(read,contigs=c("chr1","chr2"))
m <- mask_scMethrix(m, low_count = NULL, high_quantile=0.9)

cpg <- lapply(1:length(src_files)-1, function(n) {

  c <- which(rowSums(is.na(get_matrix(m)))==n)
  c <- sample(c,sample(50:100, 1))
  
})
 
m <- m[sort(unlist(cpg))]
suffix <- "_sub"

export_bed(m, path = dir, suffix=suffix)

files <- list.files (dir,full.names = TRUE)
files <- files[grepl(paste0(".*",suffix,".bedgraph$"), files,ignore.case = TRUE)]

file.rename(files, paste0(dir,"/C", 1:length(files),".bedgraph"))

files <- list.files (dir,full.names = TRUE)
files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]
files <- setdiff(files, src_files)

colData <- data.frame(phenotype = c(rep("C",2),rep("N",2)))
m <- read_beds(files=files,h5=FALSE,cov_idx=5,colData = colData)