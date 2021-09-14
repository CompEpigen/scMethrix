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
#' export_bed(scMethrix_data,path=paste0(tempdir(),"/export"))
#' @export
export_bed <- function(scm = NULL, path = NULL, suffix = NULL, verbose = TRUE, include = FALSE, na.rm = TRUE, header = FALSE) {
  
  meth <- cov <- NULL
  
  if (!is(scm, "scMethrix") || is.null(path)){
    stop("A valid scMethrix object and path needs to be supplied.", call. = FALSE)
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

# Convert \code{\link{methrix}} to \code{bsseq} object
# @details Takes \code{\link{methrix}} object and returns a \code{bsseq} object
# @param scm \code{\link{methrix}} object
# @return An object of class \code{bsseq}
# @examples
# \dontrun{
# data('methrix_data')
# methrix2bsseq(m = methrix_data)
# }
# @export
#
# export_bsseq <- function(scm) {
  # 
  # if (!is(scm, "scMethrix") || is.null(path)){
  #   stop("A valid scMethrix object and path needs to be supplied.", call. = FALSE)
  # }
  # 
  # n_samps <- nrow(SummarizedExperiment::colData(x = scm))
  # M_clean <- get_matrix(scm) * get_matrix(scm, type = "C")
  # M_clean[is.na(M_clean)] <- 0
  # assays(scm)[[2]][is.na(assays(scm)[[2]])] <- 0
  # 
  # b <- bsseq::BSseq(M = M_clean, Cov = get_matrix(scm, type = "C"), pData = colData(x = scm),
  #                   pos = rowData(x = scm)[, "start"], chr = rowData(x = scm)[, "chr"],
  #                   sampleNames = rownames(scm@colData))
#   b
# }

# Exports methrix object as bigWigs
# @param scm \code{\link{methrix}} object
# @param output_dir Output directory name where the files should be saved. Default getwd()
# @param samp_names sample names to export
# @examples
# \dontrun{
# data('methrix_data')
# write_bigwigs(m = methrix_data, output_dir = './temp')
# }
# @return NULL
# @importFrom rtracklayer export
# @export

# export_bigwigs = function(scm, output_dir = getwd(), samp_names = NULL){
  # 
  # if (!is(scm, "scMethrix") || is.null(path)){
  #   stop("A valid scMethrix object and path needs to be supplied.", call. = FALSE)
  # }
  # 
  # 
  # if (!dir.exists(output_dir)) {
  #   dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
  # }
  # 
  # mat_gr <- methrix::get_matrix(scm, type = "M", add_loci = TRUE, in_granges = TRUE)
  # 
  # seql = scm@metadata$chrom_sizes$length
  # names(seql) = scm@metadata$chrom_sizes$contig
  # 
  # all_samps = names(mcols(mat_gr))  
  # 
  # if(is.null(samp_names)){
  #   samp_names = all_samps
  # }else{
  #   samp_names = intersect(samp_names, all_samps)
  #   if(length(samp_names) == 0){
  #     stop("Incorrect sample names!")
  #   }
  # }
  # 
  # message("----------------------")
  # for(samp in samp_names){
  #   op_bw = paste0(output_dir, "/", samp, ".bw")
  #   message("*Writing ", op_bw)
  #   samp_gr = mat_gr[,samp]
  #   names(mcols(samp_gr)) = "score"
  #   samp_gr = samp_gr[!is.na(samp_gr$score)]
  #   seqlengths(samp_gr) = seql[names(seqlengths(samp_gr))]
  #   rtracklayer::export(samp_gr, con = paste0(output_dir, "/", samp, ".bw"), format="bigWig")
  # }
  # message("----------------------")
# }

# 
# export_seurat <- function(scm) {
#   
#   if (!is(scm, "scMethrix") || is.null(path)){
#     stop("A valid scMethrix object and path needs to be supplied.", call. = FALSE)
#   }
#   
#   
#   
# }
