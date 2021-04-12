# library(methylKit)
# 
# contigs <- c("chr1","chr3")
# regions <- GRanges(seqnames = "chr1", ranges = IRanges(1,10000)) 
# samples <- c("bed1","bed3")
# 
# chr <- rep(c("chr1","chr2","chr3"),each=1000)
# start <- 1:3000
# end <- 1:3000+1
# value <- rep(c(0,25,50,75,100),600)
# 
# bed1 <- data.frame(chr,start,end,value)
# bed2 <- data.frame(chr,start+1500,end+1500,value)
# names(bed2) <- c("chr","start","end","value")
# bed3 <- data.frame(chr,start+3000,end+3000,value)
# names(bed3) <- c("chr","start","end","value")
# 
# methylKit::df2tabix(bed1,"bed1_test.bed")
# methylKit::df2tabix(bed2,"bed2.bed")
# methylKit::df2tabix(bed3,"bed3.bed")