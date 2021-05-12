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
#' @param coverage flag for including a coverage matrix in the experiment
#' @param meth_idx The column index of the methylation value for the read
#' @param cov_idx The column index(es) of the read count. If not present, no coverage matrix is built
#' @export
#' @return An object of class \code{\link{scMethrix}}
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @import SingleCellExperiment GenomicRanges tools
#' @examples

# Must generate an index CpG file first:
#   sort-bed [input files] | bedops --chop 1 --ec - > CpG_index

read_beds <- function(files = NULL, ref_cpgs = NULL, colData = NULL, genome_name = "hg19", batch_size = 200, n_threads = 0, 
                      h5 = FALSE, h5_dir = NULL, h5_temp = NULL, desc = NULL, verbose = FALSE,
                      zero_based = FALSE, reads = NULL, replace = FALSE, stranded = FALSE, coverage = FALSE,
                      meth_idx = 4, cov_idx = NULL) {
  
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
  
  if (h5) {
    
    if (is.null(h5_dir)) stop("Output directory must be specified", call. = FALSE)
    
    if (is.null(ref_cpgs)) ref_cpgs <- read_index(files,n_threads,batch_size = batch_size, zero_based = zero_based)
    
    #if (zero_based) {ref_cpgs[,2:3] <- ref_cpgs[,2:3]+1}
    if (is.null(reads)) reads <- read_hdf5_data(files, ref_cpgs, n_threads, h5_temp, zero_based, verbose,
                                                meth_idx = meth_idx, cov_idx = cov_idx)
    message("Building scMethrix object")

    ref_cpgs <- GenomicRanges::makeGRangesFromDataFrame(ref_cpgs)
    colData <- data.frame()[1:(length(files)), ]
    row.names(colData) <- unlist(lapply(files,get_sample_name))
    m_obj <- create_scMethrix(methyl_mat=reads$meth, cov_mat=reads$cov, rowRanges=ref_cpgs, is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = genome_name,desc = desc,colData = colData,
                              replace = replace)
    
    message("Object built!")
    
    return(m_obj)
    
  } else {

    message("Reading in BED files") 
    
    if (is.null(ref_cpgs)) ref_cpgs <- read_index(files,n_threads,zero_based = zero_based)
    reads <- read_mem_data(files, ref_cpgs, batch_size = batch_size, n_threads = n_threads,
                           zero_based = zero_based,verbose = verbose,
                           meth_idx = meth_idx, cov_idx = cov_idx)
    ref_cpgs <- GenomicRanges::makeGRangesFromDataFrame(ref_cpgs)
    
    #colData <- t(data.frame(lapply(files,get_sample_name),check.names=FALSE))
    #colData <- t(unlist(lapply(files,get_sample_name)))
    
    message("Creating scMethrix object")
    colData <- data.frame()[1:(length(files)), ]
    row.names(colData) <- unlist(lapply(files,get_sample_name))
    m_obj <- create_scMethrix(methyl_mat=reads$meth, cov_mat=reads$cov, 
                              rowRanges=ref_cpgs, is_hdf5 = FALSE, genome_name = genome_name, 
                              desc = desc, colData = colData )
  }
}

#' Parse BED files for unique genomic coordinates
#' @details Create list of unique genomic regions from input BED files. Populates a list of batch_size+1 with 
#' the genomic coordinates from BED files, then runs \code{\link{unique}} when the list is full and keeps the running
#' results in the batch_size+1 position. Also indexes based on 'chr' and 'start' for later searching.
#' @param files List of BED files
#' @param n_threads integer for number of threads to use
#' @param batch_size Number of files to process before running unique. Default of 30.
#' @param zero_based Whether the input data is 0 or 1 based
#' @param verbose flag to output messages or not.
#' @return data.table containing all unique genomic coordinates
#' @import data.table parallel doParallel
#' @examples
#' @export
read_index <- function(files, n_threads = 0, zero_based = FALSE, batch_size = 200, verbose = TRUE) {
  
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
                                  batch_size=round(batch_size/n_threads), 
                                  n_threads = 0, zero_based = zero_based, verbose = verbose))
    
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
    data <- data.table::fread(files[i], header=FALSE, select = c(1:2))
    
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
#' @param meth_idx The column index of the methylation value for the read
#' @param cov_idx The column index(es) of the read count
#' @return data.table containing vector of all indexed methylation values for the input BED
#' @importFrom plyr .
#' @examples
read_bed_by_index <- function(file, ref_cpgs, meth_idx = 4, cov_idx = NULL, zero_based=FALSE) {
  chr <- start <- meth <- meth1 <- meth2 <- cov <- cov1 <- cov2 <- NULL
  data <- data.table::fread(file, header = FALSE, select = c(1:2,meth_idx,cov_idx))
  
  if (zero_based) data[,2] <- data[,2]+1L
  colnames(data) <- c("chr", "start", "meth",rep("cov",length(data)-3))
  data <- data.table::setkeyv(data, c("chr","start"))
  
  if (!is.null(cov_idx)) {
    
    #Collapse the coverage values
    if (length(cov_idx) > 1) {
      data[,4] <- rowSums(data[,c(cov_idx-1), with=FALSE])
    }
    data <- data[,.(chr,start,meth,cov)]
    
    #Do the search
    sample <- cbind(data[.(ref_cpgs$chr, ref_cpgs$start)][,.(meth,cov)],
                    data[.(ref_cpgs$chr, ref_cpgs$end)][,.(meth,cov)])
    colnames(sample) <- c("meth1", "cov1", "meth2","cov2")
    
    #Get the meth values from start and end CpG by weighted mean for both reads (meth*cov)
    sample[,meth1 := meth1 * cov1]
    sample[,meth2 := meth2 * cov2]
    sample[,meth := rowSums(.SD, na.rm = TRUE), .SDcols = c("meth1", "meth2")]
    sample[,cov := rowSums(.SD, na.rm = TRUE), .SDcols = c("cov1", "cov2")]
    sample[cov == 0, cov := NA] #since above line evals NA+NA as 0
    sample[,meth := meth / cov]
    sample <- sample[,.(meth,cov)]
    
    s <- nrow(ref_cpgs)-sum(is.na(sample[,(meth)]))
  } else {
    #Get the methylation value as the mean of start and end site
    sample <- rowMeans(cbind(
      data[.(ref_cpgs$chr, ref_cpgs$start)][,(meth)], 
      data[.(ref_cpgs$chr, ref_cpgs$end)][,(meth)]), na.rm=TRUE)
    s <- nrow(ref_cpgs)-sum(is.na(sample))
  }
  
  d <- nrow(data)
  
  if (s < 0.9*d) 
  {warning(paste0("Only ",round(s/d*100,1) ,"% of CpG sites in '",basename(file),"' are present in ref_cpgs"))}
  
  # Save the sample name data
  if (!is.null(cov_idx)) {
    sample <- data.table(sapply(sample, as.integer))
    sample <- list(meth = sample[,.(meth)], cov = sample[,.(cov)])
    names(sample$meth) <- get_sample_name(file)
    names(sample$cov) <- get_sample_name(file)
  } else {
    sample <- data.table(as.integer(sample))
    names(sample) <- get_sample_name(file)
    sample <- list(sample)
  }
  
  return(sample)
}

read_bed_by_index2 <- function(file,zero_based=FALSE) {
  return(read_bed_by_index(file,get("ref_cpgs", envir=globalenv()),zero_based))
}

#' Writes methylation values from input BED files into an HDF5array
#' @details Using the generated index for genomic coordinates, creates a NA-based dense matrtix of methylation
#' values for each BED file/sample. Each column contains the meth. values for a single sample.
#' @param files The BED files to parse
#' @param ref_cpgs The index of all unique coordinates from the input BED files
#' @param n_threads The number of threads to use. 0 is the default thread with no cluster built.
#' @param h5_temp The file location to store the \code{\link{RealizationSink}} object
#' @param zero_based Boolean flag for whether the input data is zero-based or not
#' @param meth_idx The column index of the methylation value for the read
#' @param cov_idx The column index(es) of the read count
#' @param verbose flag to output messages or not.
#' @return List of \code{\link{HDF5Array}}. 1 is methylation, 2 is coverage. If no cov_idx is specified, 2 will be NULL
#' @import data.table DelayedArray HDF5Array parallel doParallel
#' @examples
read_hdf5_data <- function(files, ref_cpgs, n_threads = 0, h5_temp = NULL, zero_based = FALSE, verbose = TRUE,
                           meth_idx = 4, cov_idx = NULL) {
  
  if (verbose) message("Starting HDF5 object",start_time()) 
  
  if (is.null(h5_temp)) {h5_temp <- tempdir()}
  
  # Generate the realization sinks
  dimension <- as.integer(nrow(ref_cpgs))
  colData <- as.vector(unlist(lapply(files,get_sample_name)))
  
  M_sink <- HDF5Array::HDF5RealizationSink(dim = c(dimension, length(files)),
                                           dimnames = list(NULL,colData), type = "integer",
                                           filepath = tempfile(pattern="M_sink_",tmpdir=h5_temp),
                                           name = "M", level = 6)
  
  cov_sink <- if(is.null(cov_idx)) NULL else 
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
    parallel::clusterExport(cl,list('read_bed_by_index','read_bed_by_index2','get_sample_name'))
    
    files <- split_vector(files,n_threads,by="size")
  }
  
  # Read data to the sinks
  for (i in 1:length(files)) {
    
    if (n_threads == 0) {
      if (verbose) message("   Parsing: ", get_sample_name(files[i]),appendLF=FALSE)
      
      bed <- read_bed_by_index(files[i], ref_cpgs, meth_idx = meth_idx, cov_idx = cov_idx, 
                               zero_based = zero_based)
      
      DelayedArray::write_block(block = as.matrix(bed[[1]]), viewport = grid[[i]], sink = M_sink)
      
      if (!is.null(cov_idx)) DelayedArray::write_block(block = as.matrix(bed[[2]]),
                                                       viewport = grid[[i]], sink = cov_sink)
    } else {
      if (verbose) message("   Parsing: Chunk ",i,appendLF=FALSE)
      bed <- parallel::parLapply(cl,unlist(files[i]),fun=read_bed_by_index, ref_cpgs = ref_cpgs,
                                 meth_idx = meth_idx, cov_idx = cov_idx, zero_based = zero_based)

      
      DelayedArray::write_block(block = as.matrix(cbind(lapply(bed, `[[`, 1))), 
                                viewport = grid[[i]], sink = M_sink)      
      
      
      # DelayedArray::write_block(block = as.matrix(dplyr::bind_cols(lapply(bed, `[[`, 1))), 
      #                           viewport = grid[[i]], sink = M_sink)
      if (!is.null(cov_idx)) 
        
        
        DelayedArray::write_block(block = as.matrix(cbind(lapply(data, `[[`, 2))), 
                                  viewport = grid[[i]], sink = cov_sink)
        
        # DelayedArray::write_block(block = as.matrix(dplyr::bind_cols(lapply(data, `[[`, 2))), 
        #                                             viewport = grid[[i]], sink = cov_sink)
    }
    
    rm(bed)
    if (i%%5==0) gc()
    if (verbose) message(" (",split_time(),")")
  }
  
  if (n_threads != 0) parallel::stopCluster(cl)
  if (verbose) message("Object created in ",stop_time()) 
  
  return(list(meth = M_sink, cov = cov_sink))
}

#' Writes methylation values from input BED files into an in-memory \code{\link{RangedSummarizedExperiment}}
#' @details Using the generated index for genomic coordinates, creates a NA-based dense matrtix of methylation
#' values for each BED file/sample. Each column contains the meth. values for a single sample.
#' @param files The BED files to parse
#' @param ref_cpgs The index of all unique coordinates from the input BED files
#' @param batch_size The number of files to hold in memory at once
#' @param n_threads The number of threads to use. 0 is the default thread with no cluster built.
#' @param zero_based Boolean flag for whether the input data is zero-based or not
#' @param meth_idx The column index of the methylation value for the read
#' @param cov_idx The column index(es) of the read count
#' @param verbose flag to output messages or not.
#' @return matrix of the methylation values for input BED files
#' @import parallel doParallel
#' @examples
read_mem_data <- function(files, ref_cpgs, batch_size = 200, n_threads = 0, zero_based = FALSE, 
                          meth_idx = 4, cov_idx = NULL, verbose = TRUE) {
  
  if (verbose) message("Reading BED data...",start_time()) 
  
  if (n_threads != 0) {
    # Parallel functionality
    if (verbose) message("Starting cluster with ",n_threads," threads.")
    
    cl <- parallel::makeCluster(n_threads)  
    doParallel::registerDoParallel(cl)  
    
    parallel::clusterEvalQ(cl) #, c(library(dplyr)))
    parallel::clusterExport(cl,list('read_bed_by_index','start_time','split_time','stop_time','get_sample_name'))
    
    chunk_files <- split(files, ceiling(seq_along(files)/(length(files)/n_threads)))
    
    data <- c(parallel::parLapply(cl,chunk_files,fun=read_index, 
                                  batch_size=round(batch_size/n_threads), 
                                  n_threads = 0, zero_based = zero_based, verbose = verbose))
    
    parallel::stopCluster(cl)
  
  } else {
    # Single thread functionality
    # if (verbose) message("   Parsing: Chunk ",i,appendLF=FALSE) #TODO: Get this workings
    data <- lapply(files,read_bed_by_index,ref_cpgs = ref_cpgs,zero_based = zero_based,
                   meth_idx = meth_idx, cov_idx = cov_idx)
    
    if (!is.null(cov_idx)) {
      
      data <- list(meth = cbind(lapply(data, `[[`, 1)),
                   cov = cbind(lapply(data, `[[`, 2)))
      
      # data <- list(meth = dplyr::bind_cols(lapply(data, `[[`, 1)),
      #              cov = dplyr::bind_cols(lapply(data, `[[`, 2)))
    } else {
      data <- list(meth = cbind(lapply(data, `[[`, 1)),NULL)
      
      # data <- list(meth = dplyr::bind_cols(lapply(data, `[[`, 1)),NULL)
    }
    
    # if (verbose) message(" (",split_time(),")")
  }
  
  if (verbose) message("Data read in ",stop_time()) 
  
  return (data)
}