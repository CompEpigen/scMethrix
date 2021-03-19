library(data.table)
setwd("D:\\Documents\\School\\Thesis\\scMethrix\\sample.data\\bedtools")

### Convert to tabix in bash
# From bash:
# find . -name "*.bedgraph" -exec bgzip {} \; 
# find . -name "*.gz" -exec tabix -p bed {} \;

numfiles = 1      # Number of files to generate
numrows = 1000000     # Max number of CpG sites in sample
chrs = 10          # Number of chromosomes
minsparsity = 1     # Minimum sparsity (minrows = numrows*sparsity)
maxsparsity = 1 # Max sparsity
rangeFactor = 2   # Max range of IRange (1:rangeFactor*numrows)
randomize = FALSE  # Randomize the chr and IRange mapping
values <- c(0,25,50,75,100) # List of values to choose from

start.time <- Sys.time()

  for (n in 1:numfiles) {
    
    range <- sample(x = 1:(numrows*rangeFactor), size = numrows)
    if (!randomize) range <- sort(range)
    names <- rep(c(1:chrs),each=(numrows/chrs))
    val <- sample(values, size = numrows, replace = TRUE)
    dat <- data.frame(chr = names, start = range, end = range+1, val = val)
    dat <- dat[sample(NROW(dat), NROW(dat)*(runif(1, minsparsity, maxsparsity ))),]
    dat <- dat[with(dat, order(chr, start, end)),]
    dat$chr <- paste("chr",dat$chr,sep="")
    fwrite(dat, file = paste0("bed",formatC(n, width=3, flag="0"),".bedgraph"), quote=FALSE, sep='\t', col.names = FALSE, scipen=999)
    
    message(paste("Writing:",n))

  }

message(paste0("Generating ",n," files took ",round(Sys.time() - start.time,2),"s"))
        



        