get_region_summary = function (m) {
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
}



#--------------------------------------------------------------------------------------------------------------------------
#' Loads HDF5 scMethrix object
#' @details Takes  directory with a previously saved HDF5Array format \code{\link{scMethrix}} object and loads it
#' @param dir The directory to read in from. Default NULL
#' @param ... Parameters to pass to loadHDF5SummarizedExperiment
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data('scMethrix_data')
#' methrix_data_h5 <- convert_methrix(m=methrix_data)
#' target_dir = paste0(getwd(), '/temp1/')
#' save_HDF5_methrix(methrix_data_h5, dir = target_dir, replace = TRUE)
#' load_HDF5_methrix(target_dir)
#' @export
load_HDF5_scMethrix <- function(dir = NULL, ...) {
  
  if (is.null(dir)) {
    stop("Please provide the target directory containing ")
  }
  
  m <- HDF5Array::loadHDF5SummarizedExperiment(dir = dir, ...)
  m <- as(m, "scMethrix")
  return(m)
}


#--------------------------------------------------------------------------------------------------------------------------
#' Converts HDF5 methrix object to standard in-memory object.
#' @details Takes a \code{\link{methrix}} object and returns with the same object with in-memory assay slots.
#' @param m An object of class \code{\link{methrix}}, HDF5 format
#' @return An object of class \code{\link{methrix}}
#' @examples
#' data(methrix_data)
#' m2 <- convert_methrix(m=methrix_data)
#' m <- convert_HDF5_methrix(m=m2)
#' @export
convert_HDF5_methrix <- function(m = NULL) {
  
  if (is.null(m) | !is(m, "scMethrix")) {
    stop("Invalid input data provided.")
  }
  if (!is_h5(m)) {
    stop("The input data is not in HDF5 format. Conversion aborted.")
  }

  assays(m)[[1]] <- as.matrix(assays(m)[[1]])
  m@metadata$is_h5 <- FALSE
  return(m)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Converts an in-memory object to an on-disk HDF5 object.
#' @details Takes a \code{\link{methrix}} object and returns with the same object with delayed array assay slots
#' with HDF5 backend. Might take long time!
#' @param m An object of class \code{\link{methrix}}
#' @return An object of class \code{\link{methrix}}, HDF5 format
#' @examples
#' data(methrix_data)
#' m2 <- convert_methrix(m=methrix_data)
#' @export
convert_methrix <- function(m = NULL) {
  
  if (is.null(m) | !is(m, "scMethrix")) {
    stop("No or not valid input data provided.")
  }
  if (is_h5(m)) {
    stop("The input data is already in HDF5 format. No conversion happened.")
  }
  
  m <- create_methrix(beta_mat = assays(m)[[1]],
                      cpg_loci = rowRanges(m), is_hdf5 = TRUE, genome_name = m@metadata$genome,
                      col_data = m@colData, chrom_sizes = m@metadata$chrom_sizes, ref_cpg_dt = m@metadata$ref_CpG,
                      desc = m@metadata$descriptive_stats)
  return(m)
}


  
#--------------------------------------------------------------------------------------------------------------------------
#' Order methrix object by SD
#' @details Takes \code{\link{scMethrix}} object and reorganizes the data by standard deviation
#' @param m \code{\link{scMethrix}} object
#' @param zero.rm Removes zero values from equations (the default empty value for sparse matrices)
#' @param na.rm Removes the NA values from equations
#' @return An object of class \code{\link{scMethrix}}
#' @examples
#' data('scMethrix_data')
#' order_by_sd(m = scMethrix_data)
#' @export
order_by_sd <- function (m, zero.rm = FALSE, na.rm = FALSE) {
  
#  if (!is(m, "scMethrix")){
#    stop("A valid scMethrix object needs to be supplied.")
 # }

  if (zero.rm) {
    sds <- numeric(nrow(m))
    for (i in 1:length(sds)) {
      sds[i] <- sd(as.numeric(m[i,][m[i,]!=0]), na.rm)
      print(paste(i,": ", sds[i]))
      }
  } else {
    sds <- sparseMatrixStats::rowSds(m,na.rm)
  }
  
  row_order <- order(sds, na.last = TRUE, decreasing = TRUE)
  m <- m[row_order, ]
  
  return (m)
  
}

#--------------------------------------------------------------------------------------------------------------------------

#' Subsets \code{\link{scMethrix}} object based on given conditions.
#' @details Takes \code{\link{scMethrix}} object and filters CpGs based on coverage statistics
#' @param m \code{\link{scMethrix}} object
#' @param regions genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param contigs string of chromosome names to subset by
#' @param samples string of sample names to subset by
#' @examples
#' data('scMethrix_data')
#' contigs <- c("chr1","chr3")
#' regions <- GRanges(seqnames = "chr1", ranges = IRanges(1,150)) 
#' samples <- c("bed1","bed3")
#' 
#' #Subset to chromosome 1
#' subset_scMethrix(scMethrix_data, contigs = contigs)
#' 
#' #Subset to samples bed1 and bed3
#' subset_scMethrix(scMethrix_data, samples = samples)
#' 
#' #Subset to samples bed1 and bed3, and chromosome 1
#' subset_scMethrix(scMethrix_data, samples = samples, contigs = contigs)
#' 
#' #Subset to region "chr1:1-150"
#' subset_scMethrix(scMethrix_data, regions = regions)

#' @return An object of class \code{\link{scMethrix}}
#' @export

subset_scMethrix <- function(m, regions = NULL, contigs = NULL, samples = NULL, verbose=TRUE) {
  
  message("Subsetting scMethrix")
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is_ondisk(m)) {
    
    files <- get_files(m)
    subset <- rowRanges(m)
    
    if (length(files) == 0) stop("No tabix files to subset", call. = FALSE)
    
    if (!is.null(samples)) {
      message("   Subsetting by samples")
      idx <- unlist(lapply(files,get_sample_name))
      idx <- which(idx %in% samples)   
      files <- files[idx]
      if (length(files) == 0) stop("Samples not present in files", call. = FALSE)
    }
    
    if (!is.null(contigs)) {
      message("   Subsetting by contigs")
      subset <- subset[seqnames(subset) == contigs]
      if (length(subset) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
    }
    
    if (!is.null(regions)) {
      message("   Subsetting by regions")
      regions <- cast_granges(regions)
      subset <- subsetByOverlaps(subset, regions)
      if (length(subset) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
    }
    
    subset <- sort(subset) # TODO: Somehow subset becomes unsorted; not sure why
    overlaps <- NULL
    
    for (file in files) {
      message(paste("   Parsing:",get_sample_name(file)))
      score <- get_tabix_scores(file, subset)
      if (is.null(overlaps)) {overlaps <- score
      } else {overlaps <- cbind(overlaps,score)}
    }
    
    if (nrow(overlaps) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)

    m_obj <- create_scMethrix(methyl_mat = overlaps,rowRanges = subset)
    m_obj@metadata <- m@metadata

  } else {
  
    if (!is.null(regions)) {
      message("   Subsetting by regions")
      subset <- cast_granges(regions)
      overlaps <- subsetByOverlaps(m, subset) 
      if (nrow(overlaps) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
      m_obj <- overlaps
    }
    
    if (!is.null(contigs)) {
      message("   Subsetting by contigs")
      overlaps <- m[seqnames(m) == contigs]
      if (nrow(overlaps) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
      m_obj <- overlaps
    }
    
    if (!is.null(samples)) {
      message("   Subsetting by samples")
      overlaps <- subset(m, select=colData(m)@rownames %in% samples)
      if (length(overlaps) == 0) stop("Samples not present in the object", call. = FALSE)
      m_obj <- overlaps
    }
  }
  
  return(m_obj)
  
}

get_stats <- function(m, per_chr = TRUE, verbose = TRUE) {
  
  message("Generating stats")
  
  if (!is(m, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  
  
}