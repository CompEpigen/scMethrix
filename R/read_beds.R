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
#' @import SingleCellExperiment BRGenomics GenomicRanges
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

read_beds <- function(files = NULL, colData = NULL, genome_name = "hg19", n_threads = 10, 
                      h5 = FALSE, h5_dir = NULL, desc = NULL, verbose = TRUE) {
  
  #start.time <- Sys.time()
  
  if (is.null(files)) {
    stop("Missing input files.", call. = FALSE)
  }
  
  if (h5) {

    message("Starting H5 object") 
    message("Finding unique CpG sites")  

    files_split <- split(files, ceiling(seq_along(files)/(length(files)/n_threads))) ### n_threads should be optimized

    rrng <- parallel::mclapply(files_split, generate_indexes, mc.cores = 1)#n_threads)
    
    message("Building index")
    
    rrng <- data.table::rbindlist(rrng)
    setkeyv(rrng, c("chr","start"))
    rrng <- unique(rrng)

    message("Reading data")
    
    assay <- parallel::mclapply(files_split, generate_delayed_array, index = rrng, mc.cores = 1)#n_threads)
    
    message("Creating scMethrix object")
    
    assay <- do.call(acbind, lapply(assay, acbind))

    rrng <- GenomicRanges::makeGRangesFromDataFrame(rrng)

    colData <- t(data.frame(lapply(files,get_sample_name),check.names=FALSE))

    #setHDF5DumpFile(paste0(h5_dir,"/scMethrix.h5"))
    #assay <- writeHDF5Array(assay, name="assay", chunkdim=c(nrow(assay),1), verbose=TRUE)
    
    m_obj <- create_scMethrix(methyl_mat=assay, rowRanges=rrng, is_hdf5 = TRUE, h5_dir = h5_dir,
                              genome_name = genome_name,desc = desc,colData = colData)
    
  } else {
    
    message("Reading in BED files") 
    
    if (!all(grepl("\\.(bed|bedgraph)", files))) stop("Input files must be of type bed or bedgraph.", call. = FALSE)
    
    #TODO: Replace with tabix input instead
    
    beds <- BRGenomics::import_bedGraph(files,ncores=1)
    gr <- GRangesList(beds)
    gr <- makeGRangesBRG(gr,ncores=1)
    message("Generating indexes")
    gr <- BRGenomics::mergeGRangesData(gr,ncores = 1,multiplex=TRUE)
    names(mcols(gr)) <- lapply(names(mcols(gr)),get_sample_name)
    
    rng <- c(gr, NULL, ignore.mcols=TRUE) # Remove the metadata for rowRanges input
    
    message("Creating scMethrix object")
    m_obj <- create_scMethrix(methyl_mat=mcols(gr), rowRanges=rng, is_hdf5 = h5, 
                              genome_name = genome_name, desc = desc )
    
  }
  
  #message(paste0("Reading ",length(files)," files took ",round(Sys.time() - start.time,2),"s"))
  
  return(m_obj)
  
}

#' Parse BED files for unique genomic regions
#' @details Create list of unique genomic regions from input BED files. Meant to be fed into
#' generate_delayed_array
#' @param files List of BED files
#' @return data.table containing all unique genomic regions
#' @import data.table
#' @examples
generate_indexes <- function(files) {

  rrng <- vector(mode = "list", length = length(files))
  
  for (i in 1:length(files)) {
    data <- data.table::fread(files[i], header=FALSE, select = c(1:3))
    rrng[[i]] <- data
    message(paste0("   Parsing: ",get_sample_name(files[i])))
    
  }
  
  rrng <- data.table::rbindlist(rrng)
  rrng <- unique(rrng)
  colnames(rrng) <- c("chr","start","end")

  return(rrng)
  
}

#' Parse indexed BED files into a \code{\link{DelayedArray}}
#' @details Creates a delayed array of BED files. Each column is a single sample, with known methylation
#' values places at the appropriate index from generate_indexes
#' @param files List of BED files. Column 4 must be the methylation value
#' @param index Generated index from generate_indexes. Must be generated using (at minimum) using all files
#' listed in the files parameter 
#' @return \code{\link{DelayedArray}} containing all methylation values in a dense matrix of NAs
#' @import data.table
#' @examples
generate_delayed_array <- function(files,index) {
  
  assay = NULL
  
  for (i in 1:length(files)) {
    data <- data.table::fread(files[i], header = FALSE, select = c(1:4))
    colnames(data) <- c("chr", "start", "end", "value")
    x <- index[.(data$chr, data$start), which = TRUE]
    sample <- rep(NA_integer_, nrow(index))
    sample[x] <- data[[4]]
    sample <- DelayedArray(as.matrix(sample))
    colnames(sample) <- get_sample_name(files[i])
    
    if (is.null(assay)) {assay <- sample
    } else {assay <- acbind(assay, sample)}
    
    message(paste0("   Parsing: ", get_sample_name(files[i])))
  }
  
  return(assay)
}

#' Reads in a bedgraph file and creates an indexed output for methylation values
#' @details Extracts methylation values from a sample and aligns to a generated index of methylation sites
#' @param file A BED file. Column 4 must be the methylation value
#' @param index Generated index from generate_indexes. Must be generated using (at minimum) using all files
#' listed in the files parameter 
#' @return Vector containing all methylation values and NAs for unknown sites
#' @import data.table
#' @examples
read_bdg <- function(file,index) {
  data <- data.table::fread(file, header = FALSE, select = c(1:4))
  colnames(data) <- c("chr", "start", "end", "value")
  x <- index[.(data$chr, data$start), which = TRUE]
  sample <- rep(NA_integer_, nrow(index))
  sample[x] <- data[[4]]
  colnames(sample) <- get_sample_name(file)
}


read_bed_to_hdf5 <- function(files, index, h5temp = NULL) {
  
if (is.null(h5temp)) {
  h5temp <- tempdir()
}

  dimension <- as.integer(nrow(index)/2)

  grid <- DelayedArray::RegularArrayGrid(refdim = c(dimension, length(files)),
                                         spacings = c(dimension, 1L)) 
  
  
  M_sink <- HDF5Array::HDF5RealizationSink(dim = c(dimension, length(files)),
                                           dimnames = NULL, type = "double",
                                           filepath = file.path(h5temp, paste0("M_sink_", sink_counter, ".h5")), name = "M", level = 6)
  
  
  
  
  
  
  read_chunk_to_hdf5 <- function(files, index, grid, sink, h5temp) {
    
    for (i in 1:length(files)) {
     
      bed <- read_bdg()
      
      
      DelayedArray::write_block(block = as.matrix(b$bdg[, .(beta)]),
                                viewport = grid[[i]], sink = M_sink)
      
      
       
      
      
      
      
    }
    
    
    
    
    
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


