#' Versatile BedGraph reader.
#' @details Reads BED files and generates methylation matrices.
#' Optionally arrays can be serialized as on-disk HDFS5 arrays.
#' @param files BED files containing methylation d
#' @param ref_cpgs BED files containing list of CpG sites.
#' @param stranded Default c
#' @param genome_name Name of genome. Default hg19
#' @param n_threads number of threads to use. Default 1.
#' Be-careful - there is a linear increase in memory usage with number of threads. This option is does not work with Windows OS.
#' @param on_disk Input files are of the type TABIX
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

read_beds <- function(files = NULL, colData = NULL, stranded = FALSE, genome_name = "hg19", n_threads = 1, h5 = FALSE, h5_dir = NULL, h5temp = NULL, verbose = TRUE) {
  
  #start.time <- Sys.time()
  
  if (is.null(files)) {
    stop("Missing input files.", call. = FALSE)
  }
  
  if (h5) {
    
    message("Starting H5 object") 
    
    # Create the genomic ranges
    gr <- vector(mode = "list", length = length(files))
    
    for (i in 1:length(files)) {
      
      #TODO: Parallelize fread
      data <- data.table::fread(files[i], header=FALSE, select = c(1:3))
      gr[[i]] <- data
      message(paste0("   Parsing: ",get_sample_name(files[i])))
    }
    
    gr <- data.table::rbindlist(gr)
    
    message("Generating indexes")  
    gr <- unique(gr)
    colnames(gr) <- c("chr","start","end")

    message("Writing HDF5")
    m <- as(gr, "HDF5Matrix")
    
    for (i in 1:length(files)) {
      
      data <- data.table::fread(files[i], header=FALSE, select = c(1:4))
      colnames(data) <- c("chr","start","end","value")
      x <- with(join.keys(gr, data[,1:3]), which(x %in% y))
      v <- rep(NA_integer_,nrow(gr))
      v[x] <- as.vector(unlist(data[,4]))
      v <- as(as.data.frame(v), "HDF5Matrix")
      
      m <- cbind(m, v)
      message(paste0("   Parsing: ",get_sample_name(files[i])))
    }
  
    gr <- makeGRangesFromDataFrame(gr)
    
    message("Creating scMethrix object")
    m_obj <- create_scMethrix(rowRanges=gr, files=files, on_disk = TRUE)
    
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
    m_obj <- create_scMethrix(methyl_mat=mcols(gr), rowRanges=rng, files=files, on_disk = FALSE)
    
  }
  
  #message(paste0("Reading ",length(files)," files took ",round(Sys.time() - start.time,2),"s"))
  
  return(m_obj)
  
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


