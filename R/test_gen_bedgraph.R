setwd("D:\\Documents\\School\\Thesis\\scMethrix\\test.data\\10gen")

numfiles <- 10
numrows <- 1000
chrs <- 10 # must be factor of numrows
rangeFactor <- 3 # max value is rangeFactor*numrows
  
for (n in 1:numfiles) {
  
  range <- sort(sample(x = 1:(numrows*rangeFactor), size = numrows))
  names <- rep(c(1:chrs),each=numrows/chrs)
  
  data.test <- data.frame(chr = paste("chr",names,sep=""), start = range, end = range+1, val = sample(x = 0:4)*25)
  fwrite(data.test, file=paste("bed",n,".bedgraph",sep=""), quote=FALSE, sep='\t', col.names = FALSE, scipen=999)
  
  print(n)
  
}