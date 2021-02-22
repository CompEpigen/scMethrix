#' Versatile BedGraph reader.
#' @details Reads BED files and generates methylation matrices.
#' Optionally arrays can be serialized as on-disk HDFS5 arrays.
#' @param files BED files.
#' @param stranded Default c
#' @param genome_name Name of genome. Default hg19
#' @param n_threads number of threads to use. Default 1.
#' Be-careful - there is a linear increase in memory usage with number of threads. This option is does not work with Windows OS.
#' @param on_disk Input files are of the type TABIX
#' @param verbose Be little chatty ? Default TRUE.
#' @export
#' @return An object of class \code{\link{scMethrix}}
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @import parallel
#' @import SingleCellExperiment BRGenomics
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

read_beds <- function(files = NULL, colData = NULL, stranded = FALSE, genome_name = "hg19", n_threads = 1, on_disk = NULL, verbose = TRUE) {
  
  if (on_disk) {
    
    if (!all(grepl("\\.(gz|gzip)", files))) stop("Input files must be of type gz or gzip.")
    if (!all(grepl("\\.(gz|gzip)", paste0(files,".tbi")))) stop("Input files must have corresponding *.tbi (tabix) file in the same directory")
    
    gr <- NULL
    
    for (file in files) gr <- rbind(gr, data.table::fread(file, header=FALSE,nThread=8, select = c(1:3)))
    
    gr <- unique(gr)
    colnames(gr) <- c("chr","start","end")
    
    gr <- makeGRangesFromDataFrame(gr)
    
    m_obj <- create_scMethrix(rowRanges=gr, files=files, on_disk = on_disk)
    
  } else {
  
    if (!all(grepl("\\.(bed|bedgraph)", files))) stop("Input files must be of type bed or bedgraph.")
    
    beds <- BRGenomics::import_bedGraph(files,ncores=1)
    gr <- GRangesList(beds)
    gr <- makeGRangesBRG(gr,ncores=1)
    gr <- BRGenomics::mergeGRangesData(gr,ncores = 1,multiplex=TRUE)
    
    scores <- data.frame()
    
    rng <- c(gr, NULL, ignore.mcols=TRUE)
    
    m_obj <- create_scMethrix(methyl_mat=mcols(gr), rowRanges=c(gr, NULL, ignore.mcols=TRUE), files=files, on_disk = on_disk)

  }
  
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
                      out <- rep.int(NA, length(gr))
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


