#------------------------------------------------------------------------------------------------------------
#' Exports all samples in an \code{\link{scMethrix}} objects into individual bedgraph files
#' @details The structure of the bedgraph files will be a tab-deliminated structure of:
#' Chromosome | CpG start site | CpG end site | methylation score | coverage | Additional assays (if include = TRUE)
#' 
#' If additional assays are used, and headers enabled, it is up to the user to ensure that assay names are not protected in any downstream analysis of the bedgraph files
#' @inheritParams generic_scMethrix_function
#' @param path character; the \code{\link{file.path}} of the directory to save the files
#' @param suffix character; optional suffix to add to the exported bed files 
#' @param include boolean; flag to include the values of non-standard assays in the bedgraph file
#' @param header boolean; flag to add the header onto each column
#' @param na.rm boolean; flag to remove the NA values from the output data
#' @return nothing
#' @examples
#' data('scMethrix_data')
#' export_beds(scMethrix_data,path=paste0(tempdir(),"/export"))
#' @export
export_beds <- function(scm = NULL, path = NULL, suffix = NULL, verbose = TRUE, include = FALSE, na.rm = TRUE, header = FALSE) {
  
  meth <- cov <- NULL
  
  .validateExp(scm)
  
  if (is.null(path)){
    stop("A valid path needs to be supplied.", call. = FALSE)
  }
  
  if (verbose) message("Exporting beds to ",path,start_time())
  
  dir.create(path, showWarnings = FALSE)
  
  files <- row.names(scm@colData)
  rrng <- as.data.table(rowRanges(scm))
  rrng[,c("width","strand"):=NULL]
  
  if (is.null(suffix)) suffix <- "" #TODO: Should switch to some kind of regex input
  
  for (i in 1:length(files)) {
    
    file = files[i]
    
    val <- score(scm)[, file] 
    rrng[,meth := val]
    
    if (has_cov(scm)) {
      val <- counts(scm)[, file] 
      rrng[,cov := val]
    }
    
    if (include) {
      assays <- assays(scm)
    }
    
    if (na.rm) {  out <- stats::na.omit(rrng, cols="meth", invert=FALSE)
    } else {      out <- rrng}
    
    fwrite(out, paste0(path,"/",file,suffix,".bedgraph"), append = FALSE, sep = "\t", row.names = FALSE, 
           col.names = FALSE, quote = FALSE)
    
    if (verbose) message("Exported ",i," of ",length(files)," (",split_time(), ")")
  }
  
  if (verbose) message("BEDs exported in in ",stop_time())
  
  invisible()
}

#------------------------------------------------------------------------------------------------------------
#' Converts an \code{\link{scMethrix}} object to methrix object
#' @details Removes extra slot data from an \code{\link{scMethrix}} object and changes structure to match
#' \code{\link[methrix]{methrix}} format. A 'counts' assay for coverage values must be present. 
#' Functionality not supported by methrix (e.g. reduced dimensionality)
#' will be discarded.
#' @inheritParams generic_scMethrix_function
#' @param h5_dir Location to save the methrix H5 file
#' @return a \code{\link[methrix]{methrix}} object
#' @examples
#' data('scMethrix_data')
#' # convert_to_methrix(scMethrix_data)
#' @export
export_methrix <- function(scm = NULL, h5_dir = NULL) {
  chr <- m_obj <- NULL
  
  rrng <- as.data.table(rowRanges(scm))
  rrng[,c("width","end") := NULL]
  names(rrng) <- c("chr","start","strand")
  
  chrom_size <- data.frame(contig=GenomeInfoDb::seqlevels(rowRanges(scm)),length=width(range(rowRanges(scm))))
  ref_cpgs_chr <- data.frame(chr=GenomeInfoDb::seqlevels(rowRanges(scm)),N=summary(rrng$`chr`))
  
  if (!has_cov(scm)) {
    stop("scMethrix does not contain coverage data. Cannot convert to methrix object")
  }
  
  #TODO: Need to export create_methrix function in the methrix package to use this
  if (is_h5(scm)) {
    # m_obj <- methrix::create_methrix(beta_mat = get_matrix(scm,type="score"), cov_mat = get_matrix(scm,type="counts"),
    #                                  cpg_loci = rrng[, .(chr, start, strand)], is_hdf5 = TRUE, genome_name = scm@metadata$genome,
    #                                  col_data = scm@colData, h5_dir = h5_dir, ref_cpg_dt = ref_cpgs_chr,
    #                                  chrom_sizes = chrom_sizes)#, desc = descriptive_stats)
  } else {
    # m_obj <- methrix::create_methrix(beta_mat = get_matrix(scm,type="score"), cov_mat = get_matrix(scm,type="counts"),
    #                                  cpg_loci = rrng[, .(chr, start, strand)], is_hdf5 = FALSE, 
    #                                  genome_name = scm@metadata$genome, col_data = scm@colData, 
    #                                  ref_cpg_dt = ref_cpgs_chr, chrom_sizes = chrom_sizes)#, desc = descriptive_stats)
  }
  
  return(m_obj) 
}

#' Convert \code{\link{scMethrix}} to \code{bsseq} object
#' @details Takes \code{\link{scMethrix}} object and returns a \code{bsseq} object. 
#' @param scm \code{\link{methrix}} object
#' @param m_assay matrix; the assay containing methylation scores
#' @param c_assay matrix; the assay containing count scores
#' @param path string; the path of the export directory
#' @return An object of class \code{bsseq}
#' @examples
#' \dontrun{
#' data('scMethrix_data')
#' export_bsseq(scMethrix_data)
#' }
#' @export
export_bsseq <- function(scm, m_assay = "score", c_assay="counts", path = NULL) {

  .validateExp(scm)
  
  if (!has_cov(scm)) stop("BSSeq requires a coverage matrix.", call. = FALSE)

  # if (anyNA(get_matrix(scm,m_assay)) || anyNA(get_matrix(scm,c_assay)))
  #   warning("NAs present in assay. These will be filled with zero values.")
  
  M_clean <- get_matrix(scm,m_assay) * get_matrix(scm,c_assay)
  M_clean[is.na(M_clean)] <- 0
  C_clean <- get_matrix(scm,c_assay)
  C_clean[is.na(counts(scm))] <- 0
  
  b <- bsseq::BSseq(M = M_clean, Cov = C_clean, pData = colData(scm),
                    gr = rowRanges(scm), sampleNames = rownames(colData(scm)))
  return(b)
}

#' Exports scMethrix object as bigWigs
#' @param scm \code{\link{scMethrix}} object
#' @param assay string; the assay to export. Default is "score"
#' @param path string; Output directory name where the files should be saved. Default tempdir()
#' @param samp_names string; List of sample names to export
#' @examples
#' \dontrun{
#' data('scMethrix_data')
#' export_bigwigs(scm = scMethrix_data, assay = "score", output_dir = tempdir())
#' }
#' @return NULL
#' @importFrom rtracklayer export
#' @importFrom GenomeInfoDb seqlengths
#' @export
export_bigwigs = function(scm, assay = "score", path = tempdir(), samp_names = NULL){

  .validateExp(scm)
  
  if (is.null(path)){
    stop("A valid path needs to be supplied.", call. = FALSE)
  }

  if (!dir.exists(path)) {
    dir.create(path = path, showWarnings = FALSE, recursive = TRUE)
  }
  
  message("Generating bigWig files...",start_time())
  
  mat_gr <- get_matrix(scm = scm, assay = assay, add_loci = TRUE, in_granges = TRUE)
  
  seql = width(range(rowRanges(scm)))
  names(seql) = levels(seqnames(scm))
  
  if(is.null(samp_names)){
    samp_names = names(mcols(mat_gr))
  }else{
    samp_names = intersect(samp_names, names(mcols(mat_gr)))
    if(length(samp_names) == 0) stop("Incorrect sample names! No matching samples in the experiment")
  }
  
  for(samp in samp_names){
    op_bw = paste0(path, "/", samp, ".bw")
    message("   Writing ", op_bw)
    samp_gr = mat_gr[,samp]
    names(mcols(samp_gr)) = "score"
    samp_gr = samp_gr[!is.na(samp_gr$score)]
    GenomeInfoDb::seqlengths(samp_gr) = seql[names(GenomeInfoDb::seqlengths(samp_gr))]
    rtracklayer::export(samp_gr, con = paste0(path, "/", samp, ".bw"), format="bigWig")
  }
  
  message("Files generated in ",stop_time())
}


export_seurat <- function(scm,assay="score", path = NULL) {
  
  .validateExp(scm)
  if (!has_cov(scm)) stop("Seurat requires a coverage matrix.", call. = FALSE)
  
  cnt <- counts(scm)
  rownames(cnt) <- paste0("CpG",1:nrow(cnt))
  cnt[is.na(cnt)] <- 0
  
  scr <- get_matrix(scm = scm, assay = assay)
  rownames(scr) <- paste0("CpG",1:nrow(scr))
  scr[is.na(scr)] <- 0
  
  seur <- Seurat::CreateSeuratObject(counts = cnt, meta.data = as.data.frame(colData(scm)))
  seur <- Seurat::SetAssayData(object = seur, slot = "data", new.data = scr)
  
  return(seur)
}
