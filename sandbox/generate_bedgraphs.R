#setwd("D:\\Documents\\School\\Thesis\\scMethrix\\sample.data\\100gen")

### Convert to tabix in bash
# From bash:
# find . -name "*.bedgraph" -exec bgzip {} \; find . -name "*.gz" -exec tabix -p bed {} \;

generate_bedgraph <- function (dir = NULL, numfiles = 10, numrows = 100, chrs = 10, rangeFactor = 10, randomize = FALSE) {
  
  for (n in 1:numfiles) {
    
    range <- sample(x = 1:(numrows*chrs*rangeFactor), size = chrs*numrows)
    if (!randomize) range <- sort(range)
    names <- rep(c(1:chrs),each=numrows)
    values <- sample(x = 0:4*25, size = chrs*numrows, replace = TRUE)
    
    data.test <- data.frame(chr = paste("chr",names,sep=""), start = range, end = range+1, val = values)
    fwrite(data.test, file = paste0(dir,"bed",n,".bedgraph"), quote=FALSE, sep='\t', col.names = FALSE, scipen=999)
    
    print(n)
    
  }

}