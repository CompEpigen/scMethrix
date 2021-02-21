is_ondisk = function(m) {
  return(m@metadata$is_ondisk)
}


#--------------------------------------------------------------------------------------------------------------------------
# Parse genomic regions and convert them to key'd data.table
cast_ranges <- function(regions, set.key = TRUE) {
  chr <- . <- NULL
  if (is(regions, "GRanges")) {
    target_regions <- data.table::as.data.table(x = regions)
    target_regions[, `:=`(seqnames, as.character(seqnames))]
    colnames(target_regions)[seq_len(3)] <- c("chr", "start", "end")
    if (set.key){
      data.table::setDT(x = target_regions, key = c("chr", "start", "end"))}
    target_regions <- target_regions[, .(chr, start, end)]
  } else if (is(regions, "data.frame")) {
    if (all(c("chr", "start", "end") %in% colnames(regions))){
      regions <- regions[,c("chr", "start", "end")]
    } else {
      warning("Columns with names chr, start and end are not found. Assuming that the first three columns are chr, start and end.")
    }
    target_regions <- data.table::as.data.table(x = regions)
    colnames(target_regions)[seq_len(3)] <- c("chr", "start", "end")
    target_regions <- target_regions[, .(chr, start, end)]
    target_regions[, `:=`(chr, as.character(chr))]
    target_regions[, `:=`(start, as.numeric(as.character(start)))]
    target_regions[, `:=`(end, as.numeric(as.character(end)))]
    if (set.key){
      data.table::setDT(x = target_regions, key = c("chr", "start", "end"))}
  } else {
    stop("Invalid input class for regions. Must be a data.table or GRanges object")
  }
  
  target_regions
}

cast_granges <- function(regions) {
  
  
  
  
}

#------------------------------------------------------------------------
# Gets the scores from a given tabix file for a specific region (including NA)
get_tabix_scores <- function(file, regions){
  
  if(length(tbx) = 0) stop("Input tabix file has no entries.")

  tbx = Rsamtools::scanTabix(file, param = regions)
  tbx[lengths(tbx) == 0] <- NA_character_
  tbx <- lapply(tbx,function(x) {
                  if (is.na(x)) return(NA)
                  else return (unlist(strsplit(x,"\t"))[4])
  })

  tbx2 <- do.call(rbind.data.frame, tbx)

  colnames(tbx2) <- strsplit(basename(file), "[.]")[[1]][1]
  
}

