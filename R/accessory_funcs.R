is_ondisk = function(m) {
  return(m@metadata$on_disk)
}

get_files = function(m, name_only = FALSE) {
  
  if (!is_ondisk(m)) stop("Given scMethrix object is not stored on disk.")
  
  files <- m@metadata$files
  
  if (name_only) {
    files <- unlist(lapply(files,get_sample_name))
  }
  
  return(files)
}

get_sample_name = function(s) {
  return(strsplit(basename(s), "[.]")[[1]][1])
}

chunk_granges = function(gr,factor = NA, percent = NA, num = NA) { #=NULL, percent = NULL
  
  if (length(which(is.na(c(factor,percent,num))))!=2) stop("1 argument mandatory for chunking.")
  
  if (!is.na(num)) {
    splits <- length(gr)%/%num
    splits <- num*(1:splits-1)
  }
  
  if (!is.na(percent)) factor = 100/percent
  
  if (!is.na(factor)) {
    num <- floor(length(gr)/factor)
    splits <- floor(length(gr)/num)
    splits <- 1+rep(0:(splits-1))*num
  }

  grl <- List()
  
  for (i in 1:length(splits)) {
    s <- splits[i]
    grl[[i]] <- gr[s:(s+num-1)]
  }
    
  grl[[i+1]] <- gr[(last(splits)+num):length(gr)]  
    
  #Test: all(gr == unlist(chunk_granges(gr,factor=factor,percent=percent,num=num)))

  return (GRangesList(grl))
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

#------------------------------------------------------------------------
# Replace region (text or otherwise) with Granges object
cast_granges <- function(regions) {
  if (is(regions, "GRanges")) return (regions)
  
}

#------------------------------------------------------------------------
# Gets the scores from a given tabix file for a specific region (including NA)
# Regions must be a subset of rowRanges of the associated scMethrix object
# regions <- subsetByOverlaps(rowRanges(m), regions)
get_tabix_scores <- function(file, regions=NULL){
  
  if (is.null(regions)) {tbx = Rsamtools::scanTabix(file)
  } else {tbx <- Rsamtools::scanTabix(file, param = regions)}
 
  if(length(tbx) == 0) stop(paste("Input tabix file (",file,") has no entries."), call. = FALSE)
 
  tbx[lengths(tbx) == 0] <- NA_character_
  tbx <- lapply(tbx,function(x) {
                  if (anyNA(x)) {return(NA) 
                  } else {return (unlist(strsplit(x,"\t"))[4])}
  })

  tbx <- do.call(rbind.data.frame, tbx)

  colnames(tbx) <- strsplit(basename(file), "[.]")[[1]][1]
  
  return(tbx)
}

# Not sure this actually works....
row_apply <- function(m,func,...) {
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!is_ondisk(m)) stop("Row_apply cannot be applied to non-on_disk objects.")
  
  grl <- chunk_granges(m,...)
  
  result = NULL
  
  for (gr in grl) {
    
    data <- subset_scMethrix(m,regions = gr)
    if (is.null(result)) {result <- sapply(data,func)
    } else {result <- rbind(result,sapply(data,func))}

  }
  
  return(result)  
  
}

tabix2df <- function(tbx,region=NULL) {
  
  if (is.null(regions)) {tbx = Rsamtools::scanTabix(file)
  } else {tbx <- Rsamtools::scanTabix(file, param = regions)}
  
  if(length(tbx) == 0) stop(paste("Input tabix file (",file,") has no entries."), call. = FALSE)
  
  tbx <- data.frame(lapply(tbx,function(x) {str_split_fixed(x, "\t", 4)}))
  
  colnames(tbx) <- c("chr","start","end",strsplit(basename(file), "[.]")[[1]][1])
  tbx[,2] <- as.integer(tbx[,2])
  tbx[,3] <- as.integer(tbx[,3])
  tbx[,4] <- as.integer(tbx[,4])
  
  return(df)
}