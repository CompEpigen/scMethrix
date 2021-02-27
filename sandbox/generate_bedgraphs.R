setwd("D:\\Documents\\School\\Thesis\\scMethrix\\sample.data\\100gen")

### Convert to tabix in bash
# From bash:
# find . -name "*.bedgraph" -exec bgzip {} \; 
# find . -name "*.gz" -exec tabix -p bed {} \;

numfiles = 100      # Number of files to generate
numrows = 1000000  # Max number of CpG sites in sample
chrs = 10          # Number of chromosomes
sparsity = 0.5     # Minimum sparsity (minrows = numrows*sparsity)
rangeFactor = 50   # Max range of IRange (1:rangeFactor*numrows)
randomize = FALSE  # Randomize the chr and IRange mapping

start.time <- Sys.time()

  for (n in 1:numfiles) {
    
    range <- sample(x = 1:(numrows*rangeFactor), size = numrows)
    if (!randomize) range <- sort(range)
    names <- rep(c(1:chrs),each=(numrows/chrs))
    values <- sample(x = 0:4*25, size = numrows, replace = TRUE)
    dat <- data.frame(chr = paste("chr",names,sep=""), start = range, end = range+1, val = values)
    dat <- dat[sample(NROW(dat), NROW(dat)*(runif(1, sparsity, 1))),]
    dat <- dat[with(dat, order(chr, start, end)),]
    fwrite(dat, file = paste0("bed",n,".bedgraph"), quote=FALSE, sep='\t', col.names = FALSE, scipen=999)
    
    message(paste("Writing:",n))

  }

message(paste0("Generating ",n," files took ",round(Sys.time() - start.time,2),"s"))
        
        