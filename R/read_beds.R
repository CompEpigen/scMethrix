#' Versatile BedGraph reader.
#' @details Reads BED files and generates methylation matrices.
#' Optionally arrays can be serialized as on-disk HDFS5 arrays.
#' @param files BED files containing methylation d
#' @param ref_cpgs BED files containing list of CpG sites.
#' @param stranded Default c
#' @param genome_name Name of genome. Default hg19
#' @param n_threads number of threads to use. Default 1.
#' Be-careful - there is a linear increase in memory usage with number of threads. This option is does not work with Windows OS.
#' @param h5 Should the coverage and methylation matrices be stored as 'HDF5Array'
#' @param h5_dir directory to store H5 based object
#' @param h5temp temporary directory to store hdf5
#' @param desc Description of the experiment
#' @param verbose Be little chatty ? Default TRUE.
#' @export
#' @return An object of class \code{\link{scMethrix}}
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @import SingleCellExperiment BRGenomics GenomicRanges dplyr tools parallel
#' @examples
#'\dontrun{
#'bdg_files = list.files(path = system.file('extdata', package = 'methrix'),
#'pattern = '*\\.bedGraph\\.gz$', full.names = TRUE)
#' hg19_cpgs = methrix::extract_CPGs(ref_genome = 'BSgenome.Hsapiens.UCSC.hg19')
#' meth = methrix::read_bedgraphs( files = bdg_files, ref_cpgs = hg19_cpgs,
#' chr_idx = 1, start_idx = 2, M_idx = 3, U_idx = 4,
#' stranded = FALSE, zero_based = FALSE, collapse_strands = FALSE)
#'}
#'

# Must generate an index CpG file first:
#   sort-bed [input files] | bedops --chop 1 --ec - > CpG_index

read_beds <- function(files = NULL, colData = NULL, genome_name = "hg19", n_threads = 0, 
                      h5 = FALSE, h5_dir = NULL, h5_temp = NULL, desc = NULL, verbose = FALSE,
                      zero_based = FALSE, ref_cpgs = NULL, reads = NULL, replace = FALSE) {
  
  if (is.null(files)) {
    stop("Missing input files.", call. = FALSE)
  }
  
  if (!all(grepl("\\.(bed|bedgraph)", files))) {
    stop("Input files must be of type .bed or .bedgraph", call. = FALSE)
  }
  
  if (h5) {
    
    n_threads <- min(n_threads,length(files)/2) # since cannot have multiple 1 file threads
    
    if (is.null(ref_cpgs)) ref_cpgs <- read_index(files,n_threads,zero_based = zero_based)
    
    if (zero_based) {ref_cpgs[,2:3] <- ref_cpgs[,2:3]+1}
    
    if (is.null(reads)) reads <- read_hdf5_data(files, ref_cpgs, n_threads, h5_temp, zero_based, verbose)
      
    message("Building scMethrix object")
    
    ref_cpgs <- GenomicRanges::makeGRangesFromDataFrame(ref_cpgs)
    
    colData <- t(data.frame(lapply(files,get_sample_name),check.names=FALSE))
    
    m_obj <- create_scMethrix(methyl_mat=as(reads, "HDF5Array"), rowRanges=ref_cpgs, is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = genome_name,desc = desc,colData = colData,
                              replace = replace)
    
    message("Object built!")
    
    return(m_obj)
    
  } else {
    
    message("Reading in BED files") 
    
    # beds <- lapply(files, rtracklayer::import,format = "BED")
    #  gr <- GRangesList(beds)
    #  rtracklayer::export.bed(unlist(gr),con=file.path(tempdir(),'bed1.bed'))
    
    beds <- BRGenomics::import_bedGraph(files,ncores=1)
    gr <- GRangesList(beds)
    gr <- BRGenomics::makeGRangesBRG(gr,ncores=1)
    message("Generating indexes")
    gr <- BRGenomics::mergeGRangesData(gr,ncores = 1,multiplex=TRUE)
    names(mcols(gr)) <- lapply(names(mcols(gr)),get_sample_name)
    
    rrng <- c(gr, NULL, ignore.mcols=TRUE) # Remove the metadata for rowRanges input
    
    message("Creating scMethrix object")
    m_obj <- create_scMethrix(methyl_mat=as.matrix(mcols(gr)), rowRanges=rrng, is_hdf5 = FALSE, 
                              genome_name = genome_name, desc = desc )
  }
}

#' Parse BED files for unique genomic coordinates
#' @details Create list of unique genomic regions from input BED files. Populates a list of batch_size+1 with 
#' the genomic coordinates from BED files, then runs unique() when the list is full and keeps the running
#' results in the batch_size+1 position. Also indexes based on 'chr' and 'start' for later searching.
#' @param files List of BED files
#' @param batch_size Number of files to process before running unique. Default of 30.
#' @return data.table containing all unique genomic coordinates
#' @import data.table
#' @examples
read_index <- function(files, n_threads = 0, batch_size = 200, zero_based = FALSE, verbose = TRUE) {
  
  if (n_threads != 0) {
    
    if (verbose) message("Starting cluster with ",n_threads," threads.")
    
    #no_cores <- detectCores(logical = TRUE) 
    cl <- parallel::makeCluster(n_threads)  
    registerDoParallel(cl)  
    
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
    
  if (verbose) message("Generating index",start_time())
  
  rrng <- vector(mode = "list", length = batch_size+1)
  
  for (i in 1:length(files)) {
    if (verbose) message("   Parsing: ",get_sample_name(files[i]),appendLF=FALSE)
    data <- data.table::fread(files[i], header=FALSE, select = c(1:2))
    
    if (i%%batch_size != 0) {
      rrng[[i%%batch_size]] <- data
    } else {
      rrng[[batch_size]] <- data
      data <- data.table::rbindlist(rrng)
      data <- unique(data) #data <- distinct(data)# #data <- data[!duplicated(data),]
      rrng <- vector(mode = "list", length = batch_size+1)
      rrng[[batch_size+1]] <- data
    }
    if (verbose) message(" (",split_time(),")")
  }
  
  rrng <- data.table::rbindlist(rrng)
  colnames(rrng) <- c("chr","start")
  data.table::setkeyv(rrng, c("chr","start"))
  rrng <- unique(rrng)
  
  # Remove the consecutive subsequent sites
  i = data.table::rleid(rrng$start - seq_along(rrng$start))
  i = unlist(lapply(split(i, i), seq_along))
  rrng <- rrng[i %% 2 == 1]
  
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
#' @param index The index of all unique coordinates from the input BED files
#' @return data.table containing vector of all indexed methylation values for the input BED
#' @import data.table
#' @examples
read_bed_by_index <- function(file,ref_cpgs,zero_based=FALSE) {
  data <- data.table::fread(file, header = FALSE, select = c(1:2,4))
  colnames(data) <- c("chr", "start", "value")
  if (zero_based) {data[,2] <- data[,2]+1}
  data <- data.table::setkeyv(data, c("chr","start"))
  x <- ref_cpgs[.(data$chr, data$start), which = TRUE]
  sample <- rep(NA_integer_, nrow(ref_cpgs))
  sample[x] <- data[[3]]
  sample <- as.matrix(sample)
  colnames(sample) <- get_sample_name(file)
  
  return(sample)
}

read_bed_by_index2 <- function(file,zero_based=FALSE) {
  return(read_bed_by_index(file,ref_cpgs,zero_based))
}

#' Writes methylation values from input BED files into an HDF5array
#' @details Using the generated index for genomic coordinates, creates a NA-based dense matrtix of methylation
#' values for each BED file/sample. Each column contains the meth. values for a single sample.
#' @param files The BED files to parse
#' @param ref_cpgs The index of all unique coordinates from the input BED files
#' @param n_threads The number of threads to use. 0 is the default thread with no cluster built.
#' @param zero_based
#' @param verbose
#' @param h5_temp The file location to store the RealizationSink object
#' @return HDF5Array The methylation values for input BED files
#' @import data.table DelayedArray HDF5Array
#' @examples
read_hdf5_data <- function(files, ref_cpgs, n_threads = 0, h5_temp = NULL, zero_based = FALSE, verbose = TRUE) {
  
  if (verbose) message("Starting HDF5 object",start_time()) 
  
  if (is.null(h5_temp)) {
    h5_temp <- tempdir()
  }
  
  dimension <- as.integer(nrow(ref_cpgs))
  
  M_sink <- HDF5Array::HDF5RealizationSink(dim = c(dimension, length(files)),
                                           dimnames = NULL, type = "integer",
                                           filepath = tempfile(pattern="M_sink_",tmpdir=h5_temp), name = "M", level = 6)
  
  if (n_threads == 0) {
    
    grid <- DelayedArray::RegularArrayGrid(refdim = c(dimension, length(files)),
                                           spacings = c(dimension, 1L)) 

  } else {
    
    grid <- DelayedArray::RegularArrayGrid(refdim = c(dimension, length(files)),
                                           spacings = c(dimension, n_threads)) 
    
    cl <- parallel::makeCluster(n_threads)  
    registerDoParallel(cl)  
    
    ref_cpgs <<- ref_cpgs #TODO: Figure out why this is needed and whether read_bed_by_index2 is necessary
    
    parallel::clusterEvalQ(cl, c(library(data.table)))
    parallel::clusterExport(cl,list('read_bed_by_index','read_bed_by_index2','get_sample_name',"ref_cpgs"))
    
    files <- split_vector(files,n_threads,by="size")
    
  }
    
    for (i in 1:length(files)) {
     
      if (n_threads == 0) {
        if (verbose) message("   Parsing: ", get_sample_name(files[i]),appendLF=FALSE)
        bed <- read_bed_by_index(files[i],ref_cpgs,zero_based)
        DelayedArray::write_block(block = bed, viewport = grid[[i]], sink = M_sink)s
      } else {
        if (verbose) message("   Parsing: Chunk ",i,appendLF=FALSE)
        bed <- parallel::parLapply(cl,unlist(files[i]),fun=read_bed_by_index, zero_based = zero_based)
        DelayedArray::write_block(block = dplyr::bind_cols(bed), viewport = grid[[as.integer(i)]], sink = M_sink)
      }
      
      rm(bed)
      if (i%%10==0) gc()
      if (verbose) message(" (",split_time(),")")
    }
    
  if (n_threads != 0) parallel::stopCluster(cl)

  if (verbose) message("Object created in ",stop_time()) 
  
  return(M_sink)
  
}

assignInNamespace(".multiplex_gr", ns = "BRGenomics",
                  function(data_in, field, ncores) {
                    # data must be *sorted*, base-pair resolution coverage data
                    if (all(vapply(data_in, function(x)
                      all(width(x) == 1L), logical(1L)))) {
                      data_in <- mclapply(data_in, sort, mc.cores = ncores)
                    } else {
                      warning(
                        .nicemsg(
                          "One or more inputs are not 'basepair resolution
                         GRanges' objects. Coercing them using
                         makeGRangesBRG()..."
                        ),
                        immediate. = TRUE
                      )
                      data_in <-
                        mclapply(data_in, makeGRangesBRG, mc.cores = ncores)
                    }
                    
                    # merge ranges
                    gr <- do.call(c, c(data_in, use.names = FALSE))
                    mcols(gr) <- NULL
                    gr <- unique(sort(gr))
                    
                    # (Fastest to keep these next steps separated, esp. for large datasets)
                    
                    # get dataframe of signal counts for each dataset in the output GRanges
                    idx <-
                      mclapply(data_in, function(x)
                        which(gr %in% x), mc.cores = ncores)
                    counts <- mcmapply(function(dat, idx, field) {
                      out <- rep.int(NA, length(gr))   #### <---- This line modified in this script. Reps NA instead of 0
                      out[idx] <- mcols(dat)[[field]]
                      out
                    },
                    data_in,
                    idx,
                    field,
                    mc.cores = ncores,
                    SIMPLIFY = TRUE)
                    
                    mcols(gr)[names(data_in)] <- counts
                    gr
                  })


