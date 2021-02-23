setwd("D:\\Documents\\School\\Thesis\\scMethrix\\test.data\\10gen")
library(data.table)

### Convert to tabix in bash
# From bash:
# find . -name "*.bedgraph" -exec bgzip {} \;
# find . -name "*.gz" -exec tabix -p bed {} \;

### Generate sequential bedgraphs
numfiles <- 10
numrows <- 100 #per chromosome
chrs <- 10 # must be factor of numrows
rangeFactor <- 3 # max value is rangeFactor*numrows
randomize = FALSE # set to true if genomic regions should be random
  
for (n in 1:numfiles) {
  
  range <- sample(x = 1:(numrows*chrs*rangeFactor), size = chrs*numrows)
  if (!randomize) range <- sort(range)
  names <- rep(c(1:chrs),each=numrows)
  values <- sample(x = 0:4*25, size = chrs*numrows, replace = TRUE)
  
  data.test <- data.frame(chr = paste("chr",names,sep=""), start = range, end = range+1, val = values)
  fwrite(data.test, file=paste("bed",n,".bedgraph",sep=""), quote=FALSE, sep='\t', col.names = FALSE, scipen=999)
  
  print(n)
  
}
