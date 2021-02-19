is_ondisk = function(m) {
  return(m@metadata$is_ondisk)
}


non_vect_code = function(files,verbose = FALSE,col_idx) {
  
  for (i in 1:length(files)) {
  
    if (i == 1) {
      
      data <- if(exists(data)) cbind(data,data_temp) else data_temp
      
    }
    
    
    # optimize the input window
    y <- optimal rows
    x <- optimal cols
    
    for (i in 1:(ncols(data)/y)) {  
      
      data_temp <- fread(input = files[i],
                       select = (x*(i-1)+1):min(x*i,cols), # min() catches % remainder
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE,
                       #colClasses = rep("numeric", 970),
                       col.names = col_idx,
                       data.table = TRUE)
    
      data <- data_temp %>% as.matrix() %>% dropNA() %>% cbind(data, .)
    
    }
  
  }
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