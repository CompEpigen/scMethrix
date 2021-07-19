library(profvis)

setwd("D:/Git/sampleData/Yunhee.GSE97179")
files <- list.files (getwd(),full.names = TRUE)
files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]
files <- files[1:3]
scm.big.mem <- read_beds(files=files,h5=FALSE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, M_idx=5, U_idx=6)
scm <- scm.big.mem[1:10000,]
scm.big.h5 <- read_beds(files=files,h5=TRUE,h5_dir=paste0(getwd(),"/sse"),chr_idx=1, start_idx=2, end_idx=3, 
                        beta_idx=4, M_idx=5, U_idx=6,replace=TRUE)
scm.h5 <- scm.big.h5[1:10000,]

profvis({
  #mask_by_coverage(scm,low_threshold = 1, avg_threshold = 1)
  
 # get_region_summary(scm)
  bin_scMethrix2(scm)
})