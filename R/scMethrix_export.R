#---- exportBed ------------------------------------------------------------------------------------------------------
#' Exports all samples in an [scMethrix] object into individual BED files.
#' @details The structure of the BED files will be a tab-deliminated structure of:
#' Chromosome | CpG start site | CpG end site | assay data (from those specified in `assays`) | colData (if `colData = TRUE`)
#' 
#' By default, the 4 column structure of the [BedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) format will be exported (chr-start-end-score).
#' 
#' Caution: If headers are enabled, it is up to the user to ensure that assay names or `colData()` columns are not protected in any downstream analysis of the BED files. See [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) for examples.
#' @inheritParams generic_scMethrix_function
#' @param path string; the `file.path` of the output directory. Default: `tempdir()`
#' @param suffix string; optional suffix to add to the exported bed files 
#' @param assays string; list of assays to include in each file. Each assay will be a separate column.
#' @param colData boolean; include the sample metadata from `colData()`
#' @param header boolean; flag to add the header onto each column. Some header names are invalid; see caution in details.
#' @param na.rm boolean; flag to remove the NA values from the output data
#' @param trackline boolean; flag to add track data to the file (for USCS Genome Browser)
#' @return nothing
#' @examples
#' data('scMethrix_data')
#' #exportBed(scMethrix_data,path=paste0(tempdir(),"/export"))
#' @export
exportBed <- function(scm = NULL, path = tempdir(), suffix = NULL, assays = "score", colData = FALSE, na.rm = TRUE, header = FALSE, trackline = FALSE, verbose = TRUE) {

  #---- Input validation ---------------------------------------------------
  meth <- cov <- NULL

  .validateExp(scm)
  .validateType(path,c("directory","string"))
  .validateType(suffix,c("string","null"))
  .validateType(assays,"string")
  .validateType(na.rm,"boolean")
  .validateType(header,"boolean")
  .validateType(verbose,"boolean")

  doneFiles <- NULL
  
  # Set track parameters
  parameters <-c("color" = "255,0,0",
                 "visibility"="full",
                 "altColor" = "128,128,128",
                 "autoScale"="on",
                 "viewLimits"="0:1",
                 "windowingFunction"="mean")
  parameters <- paste0(" ", paste(names(parameters), parameters,
                                  sep = "=", collapse = " "))
  
  #---- Function code ------------------------------------------------------
  if (verbose) message("Exporting beds to ",path,start_time())

  dir.create(path, showWarnings = FALSE)

  samples <- row.names(scm@colData)
  rrng <- data.table::as.data.table(rowRanges(scm))
  rrng[,c("width","strand"):=NULL]

  if (is.null(suffix)) suffix <- "" #TODO: Should switch to some kind of regex input

  for (i in 1:length(samples)) {

    samp = samples[i]

    val <- score(scm)[, samp]
    rrng[,meth := val]

    if (has_cov(scm)) {
      val <- counts(scm)[, samp]
      rrng[,cov := val]
    }

    if (include) {
      assays <- assays(scm)
    }

    if (na.rm) {  out <- stats::na.omit(rrng, cols="meth", invert=FALSE)
    } else {      out <- rrng}

    trackdata <- data.table(paste0('track type=bedGraph name="', rownames(colData(m))[i], '"', parameters))
  
    filepath <- paste0(path,"/",samp,suffix,".bedgraph")
    
    if (trackline) data.table::fwrite(trackdata, filepath, append = FALSE, sep = "\t", 
                                      row.names = FALSE,col.names = FALSE, quote = FALSE)
  
  
    data.table::fwrite(out, filepath, append = FALSE, sep = "\t", row.names = FALSE,
           col.names = FALSE, quote = FALSE)
  
    doneFiles <- c(doneFiles,filepath)
    
    if (verbose) message("Exported ",i," of ",length(samples)," (",split_time(), ")")
  }
  
  # Check to make sure all files were successfully exported
  ls <- list.files(path = path, full.names = TRUE)
  ls <- setdiff(doneFiles,ls)
  
  if (!is.empty(ls))
    stop("Error in export. Expected files are missing after export:\n",paste(ls,collapse=",\n"))
  
  if (verbose) message("BEDs exported in in ",stop_time())
}

#---- exportBedgraph --------------------------------------------------------------------------------------------------
#' Exports all samples in an [scMethrix] object into individual UCSC-compatible bedGraph files.
#' 
#' Essentially a wrapper for [exportBed()], but with mandatory format-specific options already selected.
#' @inheritParams exportBed
#' @return Nothing
#' @export
#' @examples
#' #do nothing
exportBedGraph <- function(scm = NULL, path = tempdir(), suffix = NULL, na.rm = TRUE, verbose = TRUE) {
  exportBed(scm = NULL, path = getwd(), suffix = NULL, assays = "score", colData = FALSE, na.rm = TRUE, 
             header = FALSE, trackline = TRUE, verbose = TRUE)
}


#---- exportMultiBed --------------------------------------------------------------------------------------------------
#' Exports all samples in an [scMethrix] objects into single BED file for each assay
#' @details The structure of the BED files will be a tab-deliminated structure of:
#' Chromosome | CpG start site | CpG end site | Sample Data
#' @inheritParams exportBed
#' @param assays string; the list of assays to export BED files for. Each assay will be a seperate file.
#' @param header boolean; flag to add the header onto each column. For samples, the `sampleName()` for each will be used.
#' @return nothing
#' @examples
#' data('scMethrix_data')
#' #export_beds(scMethrix_data,path=paste0(tempdir(),"/export"))
#' @export
exportMultiBed <- function(scm, path = tempdir(), suffix = NULL, assays = "score", na.rm = TRUE, header = FALSE, 
                           verbose = TRUE) {
  
  
  
  
}


#---- exportMethrix ----------------------------------------------------------------------------------------------------
#' Converts an [scMethrix] object to methrix object
#' @details Removes extra slot data from an [scMethrix] object and changes structure to match
#' [methrix::methrix] format. A `counts` assay for coverage values must be present. 
#' Functionality not supported by `methrix` (e.g. reduced dimensionality) will be discarded.
#' @inheritParams exportBed
#' @return a [methrix::methrix] object
#' @examples
#' \dontrun{#TODO: write example}
#' @export
exportMethrix <- function(scm = NULL, path = tempdir()) {

  #---- Input validation ---------------------------------------------------
  chr <- m_obj <- NULL
  
  .validateExp(scm)
  .validateType(h5_dir,"string")
  
  .validatePackageInstall("methrix")
  
  #---- Function code ------------------------------------------------------
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
    #                                  cpg_loci = rrng[, .(chr, start, strand)], is_h5 = TRUE, genome = scm@metadata$genome,
    #                                  col_data = scm@colData, h5_dir = h5_dir, ref_cpg_dt = ref_cpgs_chr,
    #                                  chrom_sizes = chrom_sizes)#, desc = descriptive_stats)
  } else {
    # m_obj <- methrix::create_methrix(beta_mat = get_matrix(scm,type="score"), cov_mat = get_matrix(scm,type="counts"),
    #                                  cpg_loci = rrng[, .(chr, start, strand)], is_h5 = FALSE, 
    #                                  genome = scm@metadata$genome, col_data = scm@colData, 
    #                                  ref_cpg_dt = ref_cpgs_chr, chrom_sizes = chrom_sizes)#, desc = descriptive_stats)
  }
  
  return(m_obj) 
}

#---- exportBSseq ------------------------------------------------------------------------------------------------------
#' Exports an [scMethrix] object as [bsseq::BSseq] object.
#' @inheritParams exportBed
#' @param scoreAssay matrix; the assay containing methylation scores
#' @param countAssay matrix; the assay containing count values
#' @return A [bsseq::BSseq] object
#' @examples
#' \dontrun{#TODO: write example}
#' @export
exportBSseq <- function(scm, scoreAssay = "score", countAssay = "counts", path = tempdir()) {

  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  if (!has_cov(scm)) stop("BSSeq requires a coverage matrix.", call. = FALSE)
  .validateAssay(scm,scoreAssay)
  .validateAssay(scm,countAssay)
  .validateType(path,"string")
  .validatePackageInstall("bsseq")
  
  # if (anyNA(get_matrix(scm,m_assay)) || anyNA(get_matrix(scm,c_assay)))
  #   warning("NAs present in assay. These will be filled with zero values.")
  
  #---- Function code ------------------------------------------------------
  M <- get_matrix(scm,scoreAssay) * get_matrix(scm,countAssay)
  M[is.na(M)] <- 0
  Cov <- get_matrix(scm,countAssay)
  Cov[is.na(counts(scm))] <- 0
  
  b <- bsseq::BSseq(M = M, Cov = Cov, pData = colData(scm),
                    gr = rowRanges(scm), sampleNames = rownames(colData(scm)))
  return(b)
}

#---- exportBigWigs ----------------------------------------------------------------------------------------------------
#' Exports an [scMethrix] object as `bigWig`.
#' @inheritParams exportBed
#' @param assay string; the assay to export. Default = "score".
#' \dontrun{#TODO: write example}
#' @return Nothing
#' @export
exportBigWigs <- function(scm, assay = "score", path = tempdir()){

  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  .validateAssay(scm,assay)
  .validateType(path,"string")
  .validatePackageInstall("rtracklayer")
  
  if (is.null(path)){
    stop("A valid path needs to be supplied.", call. = FALSE)
  }

  #---- Function code ------------------------------------------------------
  if (!dir.exists(path)) {
    dir.create(path = path, showWarnings = FALSE, recursive = TRUE)
  }
  
  message("Generating bigWig files...",start_time())
  
  mat_gr <- get_matrix(scm = scm, assay = assay, add_loci = TRUE, in_granges = TRUE)
  
  seql = GenomicRanges::width(range(SummarizedExperiment::rowRanges(scm)))
  names(seql) = levels(GenomeInfoDb::seqnames(scm))
  samples = names(mcols(mat_gr))

  for(samp in samples){
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

#---- exportSeurat -----------------------------------------------------------------------------------------------------
#' Exports an [scMethrix] object as [Seurat::Seurat]
#' @inheritParams exportBed
#' @param assay string; the assay to export. Default = "score".
#' @return A [Seurat::Seurat] object
#' @export
#' @examples
#' \dontrun{#TODO: write example}
exportSeurat <- function(scm, assay="score") {

  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  if (!has_cov(scm)) stop("Seurat requires a coverage matrix.", call. = FALSE)
  .validateAssay(scm,assay)
  .validatePackageInstall("Seurat")
    
  #---- Function code ------------------------------------------------------
  cnt <- counts(scm)
  rownames(cnt) <- paste0("CpG",1:nrow(cnt))
  cnt[is.na(cnt)] <- 0
  
  scr <- get_matrix(scm = scm, assay = assay)
  rownames(scr) <- paste0("CpG",1:nrow(scr))
  scr[is.na(scr)] <- 0
  
  seur <- Seurat::CreateSeuratObject(counts = cnt, meta.data = as.data.frame(colData(scm)))
  seur <- Seurat::SetAssayData(object = seur, slot = "data", new.data = scr)
  
  # seur <- NormalizeData(object = seur)
  # seur <- FindVariableFeatures(object = seur)
  # seur <- ScaleData(object = seur)
  # seur <- RunPCA(object = seur)
  # seur <- FindNeighbors(object = seur)
  # seur <- FindClusters(object = seur)
  # seur <- RunTSNE(object = seur)
  # DimPlot(object = seur, reduction = "tsne")
  # Set feature metadata, AKA rowData. Super intuitive, right?
  # sce.to.seurat[["RNA"]][[]] <- as.data.frame(rowData(pbmc.sce))
  
  # rownames(rowData(scm)) <- paste0("CpG",1:nrow(scm))
  
  return(seur)
}
