#' Versatile BedGraph reader.
#' @details Reads BED files and generates methylation matrices.
#' Optionally arrays can be serialized as on-disk HDFS5 arrays.
#' @param files BED files containing methylation d
#' @param ref_cpgs BED files containing list of CpG sites. Must be zero-based genome.
#' @param colData vector of input sample names
#' @param stranded Default c
#' @param genome_name Name of genome. Default hg19
#' @param batch_size Max number of files to hold in memory at once
#' @param n_threads number of threads to use. Default 1.
#' Be-careful - there is a linear increase in memory usage with number of threads. This option is does not work with Windows OS.
#' @param h5 Should the coverage and methylation matrices be stored as \code{\link{HDF5Array}}
#' @param h5_dir directory to store H5 based object
#' @param h5_temp temporary directory to store hdf5
#' @param desc Description of the experiment
#' @param verbose flag to output messages or not.
#' @param zero_based Boolean flag for whether the input data is zero-based or not
#' @param reads Manual input of reads. Typically used for testing.
#' @param replace Boolean flag for whether to delete the contents of h5_dir before saving
#' @param pipeline Default NULL. Currently supports "Bismark_cov", "MethylDackel", "MethylcTools", "BisSNP", "BSseeker2_CGmap"
#' If not known use idx arguments for manual column assignments.
#' @param chr_idx integer; column index for chromosome in bedgraph files
#' @param start_idx integer; column index for start position in bedgraph files
#' @param end_idx integer; column index for end position in bedgraph files
#' @param beta_idx integer; column index for beta values in bedgraph files
#' @param M_idx integer; column index for read counts supporting Methylation in bedgraph files
#' @param U_idx integer; column index for read counts supporting Un-methylation in bedgraph files
#' @param strand_idx integer; column index for strand information in bedgraph files
#' @param cov_idx integer; column index for total-coverage in bedgraph files
#' @export
#' @return An object of class \code{\link{scMethrix}}
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @import SingleCellExperiment GenomicRanges tools
#' @examples
#' \dontrun{
#' #Do Nothing
#' }
# Must generate an index CpG file first:
#   sort-bed [input files] | bedops --chop 1 --ec - > CpG_index

read_beds <- function(files = NULL, ref_cpgs = NULL, colData = NULL, genome_name = "hg19", batch_size = 200, n_threads = 0, 
                      h5 = FALSE, h5_dir = NULL, h5_temp = NULL, desc = NULL, verbose = TRUE,
                      zero_based = FALSE, reads = NULL, replace = FALSE, pipeline = NULL, stranded = FALSE,
                      chr_idx = NULL, start_idx = NULL, end_idx = NULL, beta_idx = NULL,
                      M_idx = NULL, U_idx = NULL, strand_idx = NULL, cov_idx = NULL) {

  if (is.null(files)) {
    stop("Missing input files.", call. = FALSE)
  }
  
  if (!all(grepl("\\.(bed|bedgraph)", files))) {
    stop("Input files must be of type .bed or .bedgraph", call. = FALSE)
  }
  
  if (n_threads > length(files)/2){
    n_threads <- min(n_threads,length(files)/2) # cannot have multiple threads with a single file being input
    warning("Too many threads specified. Each thread must have at least 2 files to process. 
            Defaulting to n_thread = ", n_threads)
  } else if (n_threads < 0) {
    n_threads <- 0
    warning("n_threads < 0. Defaulting to 0")
  }
  
  # Get the correct indexes of the input beds
  if (is.null(pipeline)) {
    message(paste0("-Preset:        Custom"))
    col_idx <- parse_source_idx(chr = chr_idx, start = start_idx, end = end_idx,
                                beta = beta_idx, cov = cov_idx, strand = strand_idx, n_meth = M_idx,
                                n_unmeth = U_idx, verbose = verbose)
  } else {
    pipeline <- match.arg(arg = pipeline, choices = c("Bismark_cov", "MethylDackel", "MethylcTools", "BisSNP", "BSseeker2_CGmap"))
    if (verbose) {
      message(paste0("-Preset:        ", pipeline))
    }
    
    col_idx <- get_source_idx(protocol = pipeline)
    
    if (any(pipeline %in% c("Bismark_cov"))) {
      if (zero_based) {
        message("*BismarkCov files are one based. You may want to re-run with zero_based=FALSE")
      }
    }
  }
  
  col_idx <<- col_idx
  #col_idx <- col_idx$col_idx
  
  if (h5) {
    
    if (is.null(h5_dir)) stop("Output directory must be specified", call. = FALSE)
    
    if(dir.exists(h5_dir) && !replace) stop("h5_dir already exists! Use 'replace=TRUE' to replace it. All 
                                            existing data in that directory will be deleted.") 
    
    if (is.null(ref_cpgs)) {
      ref_cpgs <- read_index(files = files, col_idx = col_idx, n_threads = n_threads, 
                             batch_size = batch_size, zero_based = zero_based)
    } else {
      ref_cpgs <- subset_ref_cpgs(ref_cpgs, read_index(files=files,col_idx = col_idx,n_threads,
                                                       batch_size = batch_size, zero_based = zero_based))
    }

    #if (zero_based) {ref_cpgs[,2:3] <- ref_cpgs[,2:3]+1}
    if (is.null(reads)) reads <- read_hdf5_data(files = files, ref_cpgs = ref_cpgs, col_idx = col_idx, 
                                                n_threads = n_threads, h5_temp = h5_temp, 
                                                zero_based = zero_based, verbose = verbose)
    message("Building scMethrix object")

    ref_cpgs <- GenomicRanges::makeGRangesFromDataFrame(ref_cpgs)
    
    if (is.null(colData)) colData <- data.frame()[1:(length(files)), ]
    row.names(colData) <- unlist(lapply(files,get_sample_name))
    
    chrom_size = sapply(coverage(ref_cpgs), function(x) {length(x)-x@lengths[1]})
    
    m_obj <- create_scMethrix(assays = reads, rowRanges=ref_cpgs, is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = genome_name,desc = desc,colData = colData,
                              replace = replace,chrom_size = chrom_size)
    
    message("Object built!\n")

    
  } else {

    message("Reading in BED files") 
    
    if (is.null(ref_cpgs)) ref_cpgs <- read_index(files,col_idx = col_idx, n_threads = n_threads,zero_based = zero_based)
    reads <- read_mem_data(files, ref_cpgs, col_idx = col_idx, batch_size = batch_size, n_threads = n_threads,
                           zero_based = zero_based,verbose = verbose)
    ref_cpgs <- GenomicRanges::makeGRangesFromDataFrame(ref_cpgs)
    
    #colData <- t(data.frame(lapply(files,get_sample_name),check.names=FALSE))
    #colData <- t(unlist(lapply(files,get_sample_name)))
    
    message("Building scMethrix object")
    
    if (is.null(colData)) colData <- data.frame()[1:(length(files)), ]
    row.names(colData) <- unlist(lapply(files,get_sample_name))
    
    chrom_size = sapply(coverage(ref_cpgs), function(x) {length(x)-x@lengths[1]})
    
    m_obj <- create_scMethrix(assays = reads, 
                              rowRanges=ref_cpgs, is_hdf5 = FALSE, genome_name = genome_name, 
                              desc = desc, colData = colData, chrom_size = chrom_size )
    
  }
  
  message("Object built!\n")
  return(m_obj)
}

#' Parse BED files for unique genomic coordinates
#' @details Create list of unique genomic regions from input BED files. Populates a list of batch_size+1 with 
#' the genomic coordinates from BED files, then runs \code{\link{unique}} when the list is full and keeps the running
#' results in the batch_size+1 position. Also indexes based on 'chr' and 'start' for later searching.
#' @param files List of BED files
#' @param n_threads integer for number of threads to use
#' @param batch_size Number of files to process before running unique. Default of 30.
#' @param zero_based Whether the input data is 0 or 1 based
#' @param col_idx The column indexes for the input BED files
#' @param verbose flag to output messages or not.
#' @return data.table containing all unique genomic coordinates
#' @import parallel doParallel
#' @examples
#' \dontrun{
#' #Do Nothing
#' }
#' @export
read_index <- function(files, col_idx, n_threads = 0, zero_based = FALSE, batch_size = 200, verbose = TRUE) {

  # Parallel functionality
  if (n_threads != 0) {
    
    if (n_threads > detectCores(logical = TRUE)) {
      n_threads <- detectCores(logical = TRUE)-1
      warning("Too many threads. Defaulting to n_threads =",n_threads)
    }
    
    if (verbose) message("Starting cluster with ",n_threads," threads.")
    
    cl <- parallel::makeCluster(n_threads)  
    doParallel::registerDoParallel(cl)  
    
    parallel::clusterEvalQ(cl, c(library(data.table)))
    parallel::clusterExport(cl,list('read_index','start_time','split_time','stop_time','get_sample_name'))
    
    chunk_files <- split(files, ceiling(seq_along(files)/(length(files)/n_threads)))
    
    rrng <- c(parallel::parLapply(cl,chunk_files,fun=read_index, 
                                  batch_size=round(batch_size/n_threads), col_idx = col_idx,
                                  n_threads = 0, zero_based = zero_based, verbose = FALSE))
    
    parallel::stopCluster(cl)
    
    rrng <- data.table::rbindlist(rrng)
    data.table::setkeyv(rrng, c("chr","start"))
    rrng <- unique(rrng)
    
    return(rrng)
  }
  
  # Single thread functionality  
  # Determine all unique CpG sites from input files
  if (verbose) message("Generating index",start_time())
  
  rrng <- vector(mode = "list", length = batch_size+1)

  for (i in 1:length(files)) {
    if (verbose) message("   Parsing: ",get_sample_name(files[i]),appendLF=FALSE)
    
    data <- data.table::fread(files[i], header=FALSE, select = unname(col_idx$col_idx[c("chr","start")]))
    
    # Concat the batch list if last element, otherwise save and iterate
    if (i%%batch_size == 0) {
      
      rrng[[batch_size]] <- data
      data <- data.table::rbindlist(rrng)
      data <- unique(data) #data <- distinct(data)# #data <- data[!duplicated(data),]
      rrng <- vector(mode = "list", length = batch_size+1)
      rrng[[batch_size+1]] <- data
    } else {
      rrng[[i%%batch_size]] <- data
    }
    if (verbose) message(" (",split_time(),")")
  }
  
  #Concat the final list
  rrng <- data.table::rbindlist(rrng)
  colnames(rrng) <- c("chr","start")
  data.table::setkeyv(rrng, c("chr","start"))
  rrng <- unique(rrng)
  
  # Remove the consecutive subsequent sites
  #rrng <- rrng[data.table::rowid(collapse::seqid(rrng$start)) %% 2 == 1]
  
  #Add the end position and offset if zero based
  if (zero_based) rrng$start <- rrng$start+1
  rrng$end <- rrng$start+1
  rrng <- data.table::setkeyv(rrng, c("chr","start"))
  if (verbose) message("Index generated! (",stop_time(),")")#, Avg.", (stop_time()/length(files)),")")
  
  return(rrng)
}

#' Parses BED files for methylation values using previously generated index genomic coordinates
#' @details Creates an NA-based vector populated with methlylation values from the input BED file in the
#' respective indexed genomic coordinates
#' @param file The BED file to parse
#' @param ref_cpgs The index of all unique coordinates from the input BED files
#' @param zero_based Whether the input data is 0 or 1 based
#' @param col_idx The column indexes for the input BED files
#' @return data.table containing vector of all indexed methylation values for the input BED
#' @examples
#' \dontrun{
#' #Do Nothing
#' }
#' @export
read_bed_by_index <- function(file, ref_cpgs, col_idx = NULL, zero_based=FALSE) {

  . <- chr <- start <- beta <- meth1 <- meth2 <- cov <- cov1 <- cov2 <- NULL

  col_idx <- col_idx$col_idx
  
  # Format will be: chr | start | meth | cov
  data <- data.table::fread(file, select = unname(col_idx), col.names = names(col_idx),
                                             key = c("chr", "start"))
  
  #if (zero_based) data[,start] <- data[,start]+1L
  
  cov_idx = col_idx["cov"]
  
  if (!is.na(cov_idx)) {
    
    #Do the search
    sample <- cbind(data[.(ref_cpgs$chr, ref_cpgs$start)][,.(beta,cov)],
                    data[.(ref_cpgs$chr, ref_cpgs$end)][,.(beta,cov)])
    colnames(sample) <- c("meth1", "cov1", "meth2","cov2")
    
    #Get the meth values from start and end CpG by weighted mean for both reads (meth*cov)
    sample[,meth1 := meth1 * cov1]
    sample[,meth2 := meth2 * cov2]
    sample[,beta := rowSums(.SD, na.rm = TRUE), .SDcols = c("meth1", "meth2")]
    sample[,cov := rowSums(.SD, na.rm = TRUE), .SDcols = c("cov1", "cov2")]
    sample[cov == 0, cov := NA] #since above line evals NA+NA as 0
    sample[,beta := beta / cov]
    sample <- sample[,.(beta,cov)]
    
    s <- nrow(ref_cpgs)-sum(is.na(sample[,(beta)]))
  } else {
    #Get the methylation value as the mean of start and end site
    sample <- rowMeans(cbind(
      data[.(ref_cpgs$chr, ref_cpgs$start)][,(beta)], 
      data[.(ref_cpgs$chr, ref_cpgs$end)][,(beta)]), na.rm=TRUE)
    s <- nrow(ref_cpgs)-sum(is.na(sample))
  }
  
  d <- nrow(data)
  
  if (s < 0.9*d) 
  {warning(paste0("Only ",round(s/d*100,1) ,"% of CpG sites in '",basename(file),"' are present in ref_cpgs"))}
  
  # Save the sample name data
  if (!is.null(cov_idx)) {
    sample <- data.table(sapply(sample, as.integer))
    sample <- list(meth = sample[,.(beta)], cov = sample[,.(cov)])
    names(sample$meth) <- get_sample_name(file)
    names(sample$cov) <- get_sample_name(file)
  } else {
    sample <- data.table(as.integer(sample))
    names(sample) <- get_sample_name(file)
    sample <- list(sample)
  }
  
  return(sample)
}

#' Writes values from input BED files into an in-disk \code{\link{HDF5Array}}
#' @details Using the generated index for genomic coordinates, creates a NA-based dense matrtix of methylation
#' values for each BED file/sample. Each column contains the meth. values for a single sample.
#' @param files The BED files to parse
#' @param ref_cpgs The index of all unique coordinates from the input BED files
#' @param n_threads The number of threads to use. 0 is the default thread with no cluster built.
#' @param h5_temp The file location to store the \code{\link{RealizationSink}} object
#' @param zero_based Boolean flag for whether the input data is zero-based or not
#' @param col_idx The column indexes for the input BED files
#' @param verbose flag to output messages or not.
#' @return List of \code{\link{HDF5Array}}. 1 is methylation, 2 is coverage. If no cov_idx is specified, 2 will be NULL
#' @import DelayedArray HDF5Array parallel doParallel
#' @examples
#' \dontrun{
#' #Do Nothing
#' }
read_hdf5_data <- function(files, ref_cpgs, col_idx, n_threads = 0, h5_temp = NULL, zero_based = FALSE, 
                           verbose = TRUE) {

  if (verbose) message("Starting HDF5 object",start_time()) 
  
  if (is.null(h5_temp)) {h5_temp <- tempdir()}
  
  
  has_cov <- !is.na(col_idx$col_idx["cov"])
  
  # Generate the realization sinks
  dimension <- as.integer(nrow(ref_cpgs))
  colData <- as.vector(unlist(lapply(files,get_sample_name)))
  
  M_sink <- HDF5Array::HDF5RealizationSink(dim = c(dimension, length(files)),
                                           dimnames = list(NULL,colData), type = "integer",
                                           filepath = tempfile(pattern="M_sink_",tmpdir=h5_temp),
                                           name = "M", level = 6)
  
  cov_sink <- if(!has_cov) NULL else 
              HDF5Array::HDF5RealizationSink(dim = c(dimension, length(files)),
                                             dimnames = list(NULL, colData), type = "integer",
                                             filepath = tempfile(pattern = "cov_sink_", tmpdir = h5_temp),
                                             name = "M", level = 6)
  
  # Determine the grids for the sinks
  if (n_threads == 0) {
    grid <- DelayedArray::RegularArrayGrid(refdim = c(dimension, length(files)),
                                           spacings = c(dimension, 1L)) 
  } else {
    grid <- DelayedArray::RegularArrayGrid(refdim = c(dimension, length(files)),
                                           spacings = c(dimension, n_threads)) 
    cl <- parallel::makeCluster(n_threads)  
    doParallel::registerDoParallel(cl)  
    
    parallel::clusterEvalQ(cl, c(library(data.table)))
    parallel::clusterExport(cl,list('read_bed_by_index','get_sample_name'))
    
    files <- split_vector(files,n_threads,by="size")
  }
  
  # Read data to the sinks
  for (i in 1:length(files)) {
    
    if (n_threads == 0) {
      if (verbose) message("   Parsing: ", get_sample_name(files[i]),appendLF=FALSE)
      
      bed <- read_bed_by_index(files[i], ref_cpgs, col_idx = col_idx, zero_based = zero_based)
      
      DelayedArray::write_block(block = as.matrix(bed[["meth"]]), viewport = grid[[i]], sink = M_sink)
      
      if (has_cov) DelayedArray::write_block(block = as.matrix(bed[["cov"]]),
                                                       viewport = grid[[i]], sink = cov_sink)
    } else {
      if (verbose) message("   Parsing: Chunk ",i,appendLF=FALSE)
      bed <- parallel::parLapply(cl,unlist(files[i]),fun=read_bed_by_index, ref_cpgs = ref_cpgs,
                                 col_idx = col_idx, zero_based = zero_based)

      DelayedArray::write_block(block = as.matrix(cbind(lapply(bed, `[[`, 1))), 
                                viewport = grid[[i]], sink = M_sink)      
      
      if (!is.null(col_idx$cov)) {
        DelayedArray::write_block(block = as.matrix(cbind(lapply(bed, `[[`, 2))), 
                                  viewport = grid[[i]], sink = cov_sink)
      }
    }
    
    rm(bed)
    if (i%%5==0) gc()
    if (verbose) message(" (",split_time(),")")
  }
  
  if (n_threads != 0) parallel::stopCluster(cl)
  if (verbose) message("Object created in ",stop_time()) 
  
  if (has_cov) {
    reads = list(score = M_sink, counts = cov_sink)
  } else {
    reads = list(score = M_sink)
  }
  
  reads <- lapply(reads,function(x) as(x, "HDF5Array"))
  
  return(reads)
}

#' Writes values from input BED files into an in-memory \code{\link{matrix}}
#' @details Using the generated index for genomic coordinates, creates a NA-based dense matrtix of methylation
#' values for each BED file/sample. Each column contains the meth. values for a single sample.
#' @param files The BED files to parse
#' @param ref_cpgs The index of all unique coordinates from the input BED files
#' @param batch_size The number of files to hold in memory at once
#' @param n_threads The number of threads to use. 0 is the default thread with no cluster built.
#' @param zero_based Boolean flag for whether the input data is zero-based or not
#' @param col_idx The column index for the input BED files
#' @param verbose flag to output messages or not.
#' @return matrix of the methylation values for input BED files
#' @import parallel doParallel
#' @examples
#' \dontrun{
#' #Do Nothing
#' }
read_mem_data <- function(files, ref_cpgs, col_idx, batch_size = 200, n_threads = 0, zero_based = FALSE, 
                          verbose = TRUE) {
  
  if (verbose) message("Reading BED data...",start_time()) 
  
  has_cov <- !is.na(col_idx$col_idx["cov"])
  
  if (n_threads != 0) {
    # Parallel functionality
    if (verbose) message("Starting cluster with ",n_threads," threads.")
    
    cl <- parallel::makeCluster(n_threads)  
    doParallel::registerDoParallel(cl)  
    
    parallel::clusterEvalQ(cl, c(library(data.table)))
    parallel::clusterExport(cl,list('read_bed_by_index','start_time','split_time','stop_time','get_sample_name'))
    
    chunk_files <- split(files, ceiling(seq_along(files)/(length(files)/n_threads)))
    
    reads <- c(parallel::parLapply(cl,files,fun=read_bed_by_index, ref_cpgs = ref_cpgs, 
                                   zero_based = zero_based, col_idx = col_idx))
    
    parallel::stopCluster(cl)
  
  } else {
    # Single thread functionality
    # if (verbose) message("   Parsing: Chunk ",i,appendLF=FALSE) #TODO: Get this workings
    reads <- lapply(files,read_bed_by_index,ref_cpgs = ref_cpgs,zero_based = zero_based, col_idx = col_idx)
  }
    if (has_cov) {
      reads <- list(score = cbind(lapply(reads, `[[`, 1)),
                    counts = cbind(lapply(reads, `[[`, 2)))
    } else {
      reads <- list(score = cbind(lapply(reads, `[[`, 1)))
    }
    
    reads <- lapply(reads,function(x) as.matrix(do.call(cbind, x)))
    
    # if (verbose) message(" (",split_time(),")")
  
  if (verbose) message("Data read in ",stop_time()) 
  
  return (reads)
}