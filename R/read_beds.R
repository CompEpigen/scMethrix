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
#' @import SingleCellExperiment BRGenomics GenomicRanges tictoc
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

read_beds <- function(files = NULL, colData = NULL, genome_name = "hg19", n_threads = 1, 
                                 h5 = FALSE, h5_dir = NULL, h5_temp = NULL, desc = NULL, verbose = TRUE) {

  if (is.null(files)) {
    stop("Missing input files.", call. = FALSE)
  }
  
  
  if (!all(grepl("\\.(bed|bedgraph)", files))) {
    stop("Input files must be of type .bed or .bedgraph", call. = FALSE)
  }
  
  if (h5) {

    message("Starting H5 object") 
   
    if (is.null(h5_temp)) {
      h5_temp <- tempdir()
    }
    
    index <- read_index(files)
    
    tictoc::toc()

    message("Reading data")
    
    dimension <- as.integer(nrow(index))
    
    grid <- DelayedArray::RegularArrayGrid(refdim = c(dimension, length(files)),
                                           spacings = c(dimension, 1L)) 
    
    M_sink <- HDF5Array::HDF5RealizationSink(dim = c(dimension, length(files)),
                                             dimnames = NULL, type = "integer",
                                             filepath = file.path(h5_temp, "M_sink.h5"), name = "M", level = 6)
    
    for (i in 1:length(files)) {
      cat(paste0("   Parsing: ", get_sample_name(files[i])))
      bed <- read_bed_by_index(files[i],index)
      DelayedArray::write_block(block = bed, viewport = grid[[i]], sink = M_sink)
      rm(bed)
      if (i%%10==0) gc()
      cat(paste0(" (",capture.output(tictoc::toc()),")\n"))
    }
    
    message("Data read!")
    
    message("Building scMethrix object")
    
    index <- GenomicRanges::makeGRangesFromDataFrame(index)
    
    colData <- t(data.frame(lapply(files,get_sample_name),check.names=FALSE))
    
    m_obj <- create_scMethrix(methyl_mat=as(M_sink, "HDF5Array"), rowRanges=index, is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = genome_name,desc = desc,colData = colData)
 
    file.remove(file.path(h5_temp, "M_sink.h5"))
    
    message("Object built!")
 
    tictoc::toc()
    
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
    m_obj <- create_scMethrix(methyl_mat=mcols(gr), rowRanges=rrng, is_hdf5 = FALSE, 
                              genome_name = genome_name, desc = desc )
  }
}

#' Parse BED files for unique genomic regions
#' @details Create list of unique genomic regions from input BED files. Meant to be fed into
#' generate_delayed_array
#' @param files List of BED files
#' @return data.table containing all unique genomic regions
#' @import data.table
#' @examples
read_index <- function(files) {
  
  message("Generating index")
  
  max_files <- 30 ## Test number
  
  rrng <- vector(mode = "list", length = max_files+1)
  
  for (i in 1:length(files)) {
    message(paste0("   Parsing: ",get_sample_name(files[i])))
    data <- data.table::fread(files[i], header=FALSE, select = c(1:2))
    
    if (i%%max_files != 0) {
      rrng[[i%%max_files]] <- data
    } else {
      rrng[[max_files]] <- data
      data <- data.table::rbindlist(rrng)
      data <- unique(data)
      rrng <- vector(mode = "list", length = max_files)
      rrng[[max_files+1]] <- data
    }
  }
 
  rrng <- data.table::rbindlist(rrng)
  colnames(rrng) <- c("chr","start")
  setkeyv(rrng, c("chr","start"))
  rrng <- unique(rrng)
  rrng$end <- rrng$start+1
  
  message(paste0("Index generated!"))

    return(rrng)
}

read_bed_by_index <- function(file,index) {
  data <- data.table::fread(file, header = FALSE, select = c(1:2,4))
  colnames(data) <- c("chr", "start", "value")
  x <- index[.(data$chr, data$start), which = TRUE]
  sample <- rep(NA_integer_, nrow(index))
  sample[x] <- data[[3]]
  sample <- as.matrix(sample)
  colnames(sample) <- get_sample_name(file)
  
  return(sample)
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


