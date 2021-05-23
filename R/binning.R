bin_scMethrix <- function(m, bin_size = 100000, trans = mean) {
  
  rrng <- bin_granges(rowRanges(m),bin_size = bin_size)
  
  subsetByOverlaps(rrng,rowRanges(m))
  
}