#---- exportBed ------------------------------------------------------------------------------------------------------
#' Exports all samples in an [`scMethrix`] object into individual `BED` files.
#' @details The output `BED` files will be a tab-deliminated structure of:
#' 
#' `chr` | `start` | `end` | `assay` (from those specified in `assays`) | `rowData` (if `rowData = TRUE`)
#'
#' Caution: If `colNames = TRUE`, it is up to the user to ensure that assay names or `colData()` columns are not protected in any downstream analysis of the `BED` files. See [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) for examples.
#' @inheritParams .exportBedLike
#' @inheritParams generic_scMethrix_function
#' @return nothing
#' @examples
#' data('scMethrix_data')
#' #exportBed(scMethrix_data,path=paste0(tempdir(),"/export"))
#' @export
exportBed <- function(scm = NULL, path = tempdir(), suffix = NULL, assays = "score", na.rm = TRUE, 
                      trackName = NULL, rowNames = FALSE, colNames = FALSE, rowData = FALSE, verbose = TRUE) {
  .exportBedLike(scm = scm, path = path, suffix = suffix, assays = assays, na.rm = na.rm, trackName = trackName, 
                 extension = "bed", rowNames = FALSE, colNames = FALSE, rowData = FALSE, verbose = verbose)
}

#---- exportBedgraph --------------------------------------------------------------------------------------------------
#' Exports all samples in an [`scMethrix`] object into individual UCSC-compatible `bedGraph` files.
#' 
#' @details As per [UCSC](https://genome.ucsc.edu/goldenPath/help/bedgraph.html), the output `bedGraph` files will be a tab-deliminated structure of:
#' 
#' `chr` | `start` | `end` | `score` (or other specified assay)
#' 
#' The chromosome positions should be [zero-based](https://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1).
#' @inheritParams .exportBedLike
#' @inheritParams generic_scMethrix_function
#' @return Nothing
#' @export
#' @examples
#' #do nothing
exportBedGraph <- function(scm = NULL, path = tempdir(), assay = "score", suffix = NULL, na.rm = TRUE, verbose = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  assay <- .validateAssay(scm, assay)
  if (length(assay) != 1) stop("Invalid assay. Only 1 assay can be specified.")
  
  #---- Function code ------------------------------------------------------
  .exportBedLike(scm = scm, path = path, suffix = suffix, assays = assay, rowNames = FALSE, colNames = FALSE, 
                 rowData = FALSE, na.rm = na.rm, trackName = "bedGraph", extension = "bedgraph", verbose = verbose)
}

#---- .exportBedLike ---------------------------------------------------------------------------------------------------
#' Helper function for `exportBed()` and `exportBedGraph()`
#' @inheritParams generic_scMethrix_function
#' @param path `string`; the `file.path` of the output directory. Default = `tempdir()`
#' @param suffix `string`; optional suffix to add to the exported bed files. Default = `NULL`
#' @param assays `string`; list of assays to include in each file. Each assay will be a separate column. Default = `"score"`
# @param rowData string; include the sample metadata from `colData()`
#' @param na.rm `boolean`; flag to remove the NA values from the output data. Default = `TRUE`
#' @param trackName `string`; value for `track name=` in the header line. If NULL, the header line will be excluded. Default = `NULL`
#' @param extension `string`; the extension to put on the file. Default = `"bed"`.
#' @param rowNames `boolean`; add row names to output. Default = `FALSE`.
#' @param colNames `boolean`; add column names to output. Default = `TRUE`.
#' @param rowData `boolean`; add columns for `rowData()`. Default = `FALSE`
#' @return Nothing
#' @export
#' @examples
#' data('scMethrix_data') 
.exportBedLike <- function(scm = NULL, path = tempdir(), suffix = NULL, assays = "score", na.rm = TRUE, 
                           trackName = NULL, extension = "bed", rowNames = FALSE, colNames = TRUE, rowData = FALSE, verbose = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  meth <- cov <- NULL
  
  .validateExp(scm)
  .validateType(path,c("directory","string"))
  .validateType(suffix,c("string","null"))
  .validateType(assays,"string")
  assays <- sapply(assays,function(assay) .validateAssay(scm,assay))
  .validateType(na.rm,"boolean")
  .validateType(trackName,c("string","null"))
  .validateType(rowNames,"boolean")
  .validateType(colNames,"boolean")
  .validateType(rowData,"boolean")
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
  if (verbose) {
    export <- ifelse(!is.null(trackName), paste0(trackName,"s "), "")
    message("Exporting ",export,"to ",path,start_time())
  }
  
  dir.create(path, showWarnings = FALSE)
  
  samples <- sampleNames(scm)
  mtx <- data.table::as.data.table(granges(rowRanges(scm), use.mcols=rowData))
  mtx[,c("width","strand") := NULL]
  setnames(mtx, "seqnames", "chr")
  
  if (is.null(suffix)) suffix <- "" #TODO: Should switch to some kind of regex input
  
  for (i in 1:length(samples)) {
    
    out <- mtx
    
    samp = samples[i]
    filepath <- paste0(path,"/",samp,suffix,".",extension)
    
    for (assay in assays) {
      out[,(assay) := get_matrix(scm, assay)[, samp]]
    }

    if (na.rm) out <- out[complete.cases(out), , ]
    
    appnd <- FALSE
    
    if (!is.null(trackName)) {
      trackName <- data.table(paste0('track type=',trackName,' name="', sampleNames(scm)[i], '"', parameters))
      data.table::fwrite(trackName, filepath, append = FALSE, sep = "\t", 
                                      row.names = FALSE, col.names = FALSE, quote = FALSE)
      appnd <- TRUE
    }
    
    data.table::fwrite(out, filepath, append = appnd, sep = "\t", row.names = rowNames,
                       col.names = colNames, quote = FALSE, na = NA)
    
    doneFiles <- c(doneFiles,filepath)
    
    if (verbose) message("Exported sample ",get_sample_name(samp), " (",i," of ",length(samples),"; ",split_time(), ")")
  }
  
  # Check to make sure all files were successfully exported
  ls <- list.files(path = path, full.names = TRUE)
  ls <- setdiff(doneFiles,ls)
  
  if (!is.empty(ls))
    stop("Error in export. Expected files are missing after export:\n",paste(ls,collapse=",\n"))
  
  if (verbose) message("BEDs exported in in ",stop_time())

}

#---- exportMultiBed --------------------------------------------------------------------------------------------------
#' Exports all samples in an [`scMethrix`] objects into single `BED` file for each assay
#' @details The structure of the `MultiBED` files will be a tab-deliminated structure of:
#' `chr` | `start` | `end` | `assay`
#' @inheritParams .exportBedLike
#' @param assays `string`; the list of assays to export `BED` files for. Each assay will be a separate file.
#' @return nothing
#' @examples
#' data('scMethrix_data')
#' #export_beds(scMethrix_data,path=paste0(tempdir(),"/export"))
#' @export
exportMultiBed <- function(scm, path = tempdir(), suffix = NULL, assays = "score", na.rm = TRUE, rowNames = FALSE, colNames = TRUE, verbose = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  meth <- cov <- NULL
  
  .validateExp(scm)
  .validateType(path,c("directory","string"))
  .validateType(suffix,c("string","null"))
  assays <- sapply(assays,function(assay) .validateAssay(scm,assay))
  .validateType(assays,"string")
  .validateType(na.rm,"boolean")
  .validateType(verbose,"boolean")
  .validateType(rowNames,"boolean")
  .validateType(colNames,"boolean")
  
  doneFiles <- NULL

  #---- Function code ------------------------------------------------------
  if (verbose) message("Exporting beds to ",path,start_time())
  
  dir.create(path, showWarnings = FALSE)
  
  if (is.null(suffix)) suffix <- "" #TODO: Should switch to some kind of regex input

  assayNames <- assays # to avoid confusion with assays()
  
  for (i in 1:length(assayNames)) {

    mtx <- get_matrix(scm,assayNames[i])
    mtx <- as.data.table(mtx)
    name <- names(assays(scm))[[i]]
    
    filepath <- paste0(path,"/",name,suffix,".bed")
    
    if (na.rm) mtx <- mtx[complete.cases(mtx), , ]
    
    data.table::fwrite(mtx, filepath, append = FALSE, sep = "\t", row.names = rowNames,
                       col.names = colNames, quote = FALSE, na = NA)
    
    if (verbose) message("Exported assay '",assayNames[i], "' (",i," of ",length(assayNames),"; ",split_time(), ")")
  }
}


#---- exportMethrix ----------------------------------------------------------------------------------------------------
#' Converts an [`scMethrix`] object to a [`methrix`][methrix::methrix-class] object
#' @details Removes extra slot data from an [`scMethrix`] object and changes structure to match
#' [`methrix`][methrix::methrix-class] format. A `counts` assay for coverage values must be present. 
#' Functionality not supported by [`methrix`][methrix::methrix-class] (e.g. reduced dimensionality) will be discarded.
#' @inheritParams exportBed
#' @return a [`methrix::methrix-class`] object
#' @examples
#' \dontrun{#TODO: write example}
#' @export
exportMethrix <- function(scm = NULL, path = tempdir()) {

  #---- Input validation ---------------------------------------------------
  chr <- m_obj <- NULL
  
  .validateExp(scm)
  .validateType(path,"string")
  
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
#' Exports an [`scMethrix-class`] object as [`bsseq::BSseq`] object.
#' @inheritParams exportBed
#' @param scoreAssay `matrix`; the assay containing methylation scores
#' @param countAssay `matrix`; the assay containing count values
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
  
   if (anyNA(get_matrix(scm,scoreAssay)) || anyNA(get_matrix(scm,countAssay)))
     warning("NAs present in assay. These will be filled with zero values.")
  
  #---- Function code ------------------------------------------------------
  M <- get_matrix(scm,scoreAssay) * get_matrix(scm,countAssay)
  M[is.na(M)] <- 0
  Cov <- get_matrix(scm,countAssay)
  Cov[is.na(counts(scm))] <- 0
  
  b <- bsseq::BSseq(M = M, Cov = Cov, pData = colData(scm),
                    gr = rowRanges(scm), sampleNames = rownames(colData(scm)))
  return(b)
}

#---- exportBigWig ----------------------------------------------------------------------------------------------------
#' Exports an [`scMethrix`] object as `bigWig`.
#' @inheritParams exportBed
#' @inheritParams generic_scMethrix_function
#' @return Nothing
#' @examples 
#' \dontrun{#TODO: write example}
#' @export
exportBigWig <- function(scm, assay = "score", path = tempdir()){

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
#' Exports an [`scMethrix`] object as [Seurat::Seurat]
#' @inheritParams .exportBedLike
#' @inheritParams generic_scMethrix_function
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
  .validatePackageInstall("SeuratObject")
    
  #---- Function code ------------------------------------------------------
  cnt <- counts(scm)
  rownames(cnt) <- paste0("CpG",1:nrow(cnt))
  cnt[is.na(cnt)] <- 0
  
  scr <- get_matrix(scm = scm, assay = assay)
  rownames(scr) <- paste0("CpG",1:nrow(scr))
  scr[is.na(scr)] <- 0
  
  seur <- Seurat::CreateSeuratObject(counts = cnt, meta.data = as.data.frame(colData(scm)))
  seur <- Seurat::SetAssayData(object = seur, slot = "data", new.data = scr)
  seur <- SeuratObject::RenameAssays(seur, RNA = "Methylation")
  
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
