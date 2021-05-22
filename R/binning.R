bin_scMethrix <- function(m, bin_size = 100000) {
  
  ends <- len <- seqnames(m)@lengths
  for (i in 1:length(ends)) ends[i] <- sum(as.vector(len[1:i]))
  starts <- head(c(1, ends + 1), -1)
  
  rngs <- lapply(1:length(starts), function (i) {
    
    #Get span of each chr     
    gr <- c(rowRanges(m)[starts[i]],rowRanges(m)[ends[i]])
    gr <- c(gr,gaps(gr,start=rowRanges(m)[starts[i]]@ranges@start))
    gr <- reduce(gr)
    
    #Create bins
    tiles <- tile(gr, width = bin_size)
    
    
    
    
    
  })
  
  
  
  
  
  rowRanges(m)[starts[2]]
  rowRanges(m)[ends[2]]
  
  
  
  
  
}