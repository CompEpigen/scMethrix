#' Versatile BedGraph reader.
#' @details Reads BED files and generates methylation matrices.
#' Optionally arrays can be serialized as on-disk HDFS5 arrays.
#' @param files BED files.
#' @param stranded Default c
#' @param genome_name Name of genome. Default hg19
#' @param n_threads number of threads to use. Default 1.
#' Be-careful - there is a linear increase in memory usage with number of threads. This option is does not work with Windows OS.
#' @param h5 Should the coverage and methylation matrices be stored as 'HDF5Array'
#' @param h5_dir directory to store H5 based object
#' @param h5temp temporary directory to store hdf5
#' @param verbose Be little chatty ? Default TRUE.
#' @export
#' @return An object of class \code{\link{scMethrix}}
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @import parallel
#' @import DelayedMatrixStats
#' @import SingleCellExperiment DelayedArray HDF5Array
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

read_beds <- function(files = NULL, stranded = FALSE, genome_name = "hg19", n_threads = 1, h5 = NULL, h5_dir = NULL, h5temp = NULL, verbose = TRUE) {

  genome <- Seqinfo(genome = NA_character_)
  i <- lapply(files,basename)
  i <- file_path_sans_ext(i)
  names(i) <- samples
  bl <- List(lapply(i, import, genome = genome_name))
  gr_b <- stack(bl, "i")
  dj <- disjoin(gr_b, ignore.strand = !stranded, with.revmap = TRUE)
  mcols(gr_b)$i <- decode(mcols(gr_b)$i)
  dfl <- extractList(mcols(gr_b), mcols(dj)$revmap)
  assay <- as.matrix(dfl[, "score"], col.names = dfl[,"i"])
  rowData <- granges(dj)
  ans <- SummarizedExperiment(list(score = assay), rowData)
  sse <- as(ans, "SingleCellExperiment")





  return(sse)

}



