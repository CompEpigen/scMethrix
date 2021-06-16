#------------------------------------------------------------------------------------------------------------
#' Adds descriptive statistics to metadata columns in an \code{\link{scMethrix}} object.
#' @details Adds the mean, median and SD for each region in an \code{\link{scMethrix}} object
#' @inheritParams generic_scMethrix_function
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' get_metadata_stats(scMethrix_data)
#' @export
get_metadata_stats <- function(scm) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is_h5(scm)) {
    stats <-
      data.table::data.table(
        mean_meth = DelayedMatrixStats::rowMeans2(get_matrix(scm), na.rm = TRUE),
        median_meth = DelayedMatrixStats::rowMedians(get_matrix(scm), na.rm = TRUE),
        sd_meth = DelayedMatrixStats::rowSds(get_matrix(scm), na.rm = TRUE)
      )
    
    if(has_cov(scm)) stats[,"counts" := DelayedMatrixStats::rowSums2(get_matrix(scm,assay="counts"), na.rm = TRUE)] 
    
  } else {
    stats <- data.table::data.table(
        mean_meth = matrixStats::rowMeans2(get_matrix(scm), na.rm = TRUE),
        median_meth = matrixStats::rowMedians(get_matrix(scm), na.rm = TRUE),
        sd_meth = matrixStats::rowSds(get_matrix(scm), na.rm = TRUE)
    )
    
    if(has_cov(scm)) stats[,"counts" := matrixStats::rowSums2(get_matrix(scm,assay="counts"), na.rm = TRUE)] 
  }
  mcols(scm) <- stats
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Removes an assay from an \code{\link{scMethrix}} object
#' @inheritParams generic_scMethrix_function
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' remove_assay(scMethrix_data,assay="counts")
#' @export
remove_assay <- function(scm,assay) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay is not in object.", call. = FALSE)
  }
  
  if (assay == "score") {
    stop("Score assay cannot be removed.", call. = FALSE)
  }
  
  assays(scm) <- assays(scm)[-which(SummarizedExperiment::assayNames(scm) == assay)]
  
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Merges two \code{\link{scMethrix}} objects by \code{row} or \code{col}
#' @details Merges the base assay data from two \code{\link{scMethrix}} objects. Merging of additional slot
#' data is not supported at this time. Non-common assays between objects will be dropped
#' @param scm1 A \code{\link{scMethrix}} object
#' @param scm2 A \code{\link{scMethrix}} object
#' @param by Merge by columns or rows
#' @return A merged \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' merge_scMethrix(scMethrix_data[1:5],scMethrix_data[6:10],by="row")
#' @export
merge_scMethrix <- function(scm1 = NULL, scm2 = NULL, by = c("row", "col")) {
  
  if (!is(scm1, "scMethrix") || !is(scm2, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  names1 = SummarizedExperiment::assayNames(scm1)
  names2 = SummarizedExperiment::assayNames(scm2)
  
  if (!all((sort(names1)==sort(names2)))) {
    warning("Assay list not identical. All non-identical assays will be dropped from merged object.")
    a1 <- intersect(names1, names2)
    a2 <- intersect(names2, names1)
    assays(scm1) <- assays(scm1)[a1]
    assays(scm2) <- assays(scm2)[a2]
  } 
  
  by <- match.arg(arg = by, choices = c("row", "col"), several.ok = FALSE)
  
  if (by == "row") {
    if (nrow(colData(scm1)) != nrow(colData(scm2)) || !all(rownames(scm1@colData) == rownames(scm2@colData))) {
      stop("You have different samples in your dataset. You need the same samples in your datasets. ")
    } else {
      m <- rbind(scm1, scm2)
    }
    if (any(duplicated(rowRanges(m)))) {
      stop("There are overlapping regions in your datasets. Each object must contain unique regions. ")
    }
  }
  if (by == "col") {
    if (any(rownames(scm1@colData) %in% rownames(scm2@colData))) {
      stop("You have the same samples in your datasets. You need different samples for this merging.  ")
    } else if (!identical(rowRanges(scm1),rowRanges(scm2))) {
      stop("There are non-overlapping regions in your datasets. This function only takes identical regions. ")
    } else {
      m <- cbind(scm1, scm2)
    }
  }
  
  return(m)
}

#------------------------------------------------------------------------------------------------------------
#' Converts an \code{\link{scMethrix}} object to methrix object
#' @details Removes extra slot data from an \code{\link{scMethrix}} object and changes structure to match
#' \code{\link[methrix]{methrix}} format
#' @inheritParams generic_scMethrix_function
#' @param h5_dir Location to save the methrix H5 file
#' @return a \code{\link[methrix]{methrix}} object
#' @examples
#' data('scMethrix_data')
#' # convert_to_methrix(scMethrix_data)
#' @export
convert_to_methrix <- function(scm = NULL, h5_dir = NULL) {
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

#------------------------------------------------------------------------------------------------------------
#' Exports all samples in an \code{\link{scMethrix}} objects into individual bedgraph files
#' @inheritParams generic_scMethrix_function
#' @param path the \code{\link{file.path}} of the directory to save the files
#' @param suffix optional suffix to add to the exported bed files 
#' @return nothing
#' @examples
#' data('scMethrix_data')
#' export_bed(scMethrix_data,path=paste0(tempdir(),"/export"))
#' @export
export_bed <- function(scm = NULL, path = NULL, suffix = NULL) {
  meth <- cov <- NULL
  if (!is(scm, "scMethrix") || is.null(path)){
    stop("A valid scMethrix object and path needs to be supplied.", call. = FALSE)
  }
  
  dir.create(path, showWarnings = FALSE)
  
  files <- row.names(scm@colData)
  rrng <- as.data.table(rowRanges(scm))
  rrng[,c("width","strand"):=NULL]
  
  if (is.null(suffix)) suffix <- "" #TODO: Should switch to some kind of regex input
  
  for (file in files) {
    
    val <- as.data.table(get_matrix(scm,assay="score"))[, file, with=FALSE] 
    rrng[,meth := val]
    
    if (has_cov(scm)) {
      val <- as.data.table(get_matrix(scm,assay="counts"))[, file, with=FALSE] 
      rrng[,cov := val]
    }
    
    write.table(rrng[which(!is.na(val)),], paste0(path,"/",file,suffix,".bedgraph"), append = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  invisible()
}

#------------------------------------------------------------------------------------------------------------
#' Extracts and summarizes methylation or coverage info by regions of interest
#' @details Takes \code{\link{scMethrix}} object and summarizes regions
#' @inheritParams generic_scMethrix_function
#' @param regions genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GenomicRanges}} object
#' @param type matrix which needs to be summarized. Could be `M`, `C`. Default 'M'
#' @param how mathematical function by which regions should be summarized. Can be one of the following: mean, sum, max, min. Default 'mean'
#' @param overlap_type defines the type of the overlap of the CpG sites with the target region. Default value is `within`. For detailed description,
#' see the \code{findOverlaps} function of the \code{\link{IRanges}} package.
#' #param elementMetadata.col columns in \code{\link{scMethrix}}@elementMetadata which needs to be summarised. Default = NULL.
#' @param n_chunks Number of chunks to split the \code{\link{scMethrix}} object in case it is very large. Default = 1.
#' @param n_threads Number of parallel instances. \code{n_cores} should be less than or equal to \code{n_chunks}. If \code{n_chunks} is not specified, then \code{n_chunks} is initialized to be equal to \code{n_cores}. Default = 1.
#' @param group a column name from sample annotation that defines groups. In this case, the number of min_samples will be
#' tested group-wise.
#' @return table of summary statistic for the given region
#' @examples
#' data('scMethrix_data')
#' get_region_summary(scMethrix_data,
#'    regions = data.table(chr = c('chr1','chr2'), start = c(1), end =  c(10000000)),
#'    type = 'score', how = 'mean')
#' @export
get_region_summary = function (scm = NULL, regions = NULL, n_chunks=1, n_threads = 1, type="score", how = "mean", 
                               overlap_type = "within", verbose = TRUE, group = NULL) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!is.null(group) && !(group %in% colnames(scm@colData))){
    stop(paste("The column name ", group, " can't be found in colData. Please provid a valid group column."))
  }
  
  if (n_chunks > nrow(scm)) {
    n_chucks <- nrow(scm)
    warning("n_chunks exceeds number of files. Defaulting to n_chunks = ",n_chunks)
  }
  
  if(verbose) message("Generating region summary...",start_time())
  
  type = match.arg(arg = type, choices = SummarizedExperiment::assayNames(scm))
  how = match.arg(arg = how, choices = c('mean', 'median', 'max', 'min', 'sum', 'sd'))
  yid  <- NULL
  
  regions = cast_granges(regions)
  regions$rid <- paste0("rid_", 1:length(regions))
  
  overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, regions, type = overlap_type)) #GenomicRanges::findOverlaps(rowRanges(m), regions)@from
  
  if(nrow(overlap_indices) == 0){
    stop("No overlaps detected")
  }
  
  if(nrow(overlap_indices) == 0){
    stop("No overlaps detected")
  }
  
  colnames(overlap_indices) <- c("xid", "yid")
  overlap_indices[,yid := paste0("rid_", yid)]
  n_overlap_cpgs = overlap_indices[, .N, yid]
  colnames(n_overlap_cpgs) = c('rid', 'n_overlap_CpGs')
  
  if(n_chunks==1){
    if (type == "score") {
      dat = get_matrix(scm = scm[overlap_indices$xid,], assay = "score", add_loci = TRUE)
    } else if (type == "counts") {
      dat = get_matrix(scm = scm[overlap_indices$xid,], assay = "counts", add_loci = TRUE)
    }
  } else {
    
    stop("Chunking not enabled")
    
    # if(nrow(overlap_indices) < n_chunks){
    #   n_chunks <- nrow(overlap_indices)
    #   warning("Fewer overlaps indicies than n_chunks. Defaulting to n_chunks = ",n_chunks)
    # }
    # 
    # if (n_chunks < n_threads) {
    #   n_threads <- n_chunks
    #   warning("n_threads < n_chunks. Defaulting to n_threads = ",n_threads)
    # }
    # 
    # cl <- parallel::makeCluster(n_threads)  
    # doParallel::registerDoParallel(cl)  
    # 
    # parallel::clusterEvalQ(cl, c(library(data.table), sink(paste0("D:/Git/scMethrix/", Sys.getpid(), ".txt"))))
    # parallel::clusterExport(cl,list('m','scMethrix','type', 'get_matrix','start_time','split_time','stop_time'))
    # 
    # chunk_overlaps <- split(overlap_indices$xid, ceiling(seq_along(overlap_indices$xid) /
    #                                                     ceiling(length(overlap_indices$xid)/n_chunks)))
    # 
    # data = c(parallel::parLapply(cl,chunk_overlaps,fun=function(i) {
    #   cat("Looking for i =",i,typeof(i))
    #   get_matrix(m[i,], type = type, add_loci = TRUE) # TODO: object of type 'S4' is not subsettable
    # }))
    # 
    # # data = lapply(chunk_overlaps,FUN=function(i) {
    # #   cat("Looking for i =",toString(i))
    # #   get_matrix(m[i,], type = type, add_loci = TRUE)
    # # })
    # 
    # parallel::stopCluster(cl)
    # 
    # data <- rbindlist(data)
    
  }
  
  if(nrow(overlap_indices) != nrow(dat)){
    stop("Something went wrong")
  }
  
  dat = cbind(overlap_indices, dat)
  
  #message("-Summarizing overlaps..\n")
  if(how == "mean") {
    message("Summarizing by average")
    
    output = dat[, lapply(.SD, mean, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(scm)))]
  } else if (how == "median") {
    message("Summarizing by median")
    output = dat[, lapply(.SD, median, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(scm)))]
  } else if (how == "max") {
    message("Summarizing by maximum")
    output = suppressWarnings(
      dat[, lapply(.SD, max, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(scm)))])
    
    for (j in 1:ncol(output)) set(output, which(is.infinite(output[[j]])), j, NA)
  } else if (how == "min") {
    message("Summarizing by minimum")
    output = suppressWarnings(
      
      dat[, lapply(.SD, min, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(scm)))])
    for (j in 1:ncol(output)) set(output, which(is.infinite(output[[j]])), j, NA)
  } else if (how == "sum") {
    message("Summarizing by sum")
    output = dat[, lapply(.SD, sum, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(scm)))]
  }
  
  output = merge(regions, output, by.x = 'rid', by.y = 'yid', all.x = TRUE)
  
  output = merge(n_overlap_cpgs, output, by = 'rid')
  output$rid <- as.numeric(gsub("rid_","",output$rid))
  
  output <- output[order(output$rid),]
  setnames(output, "seqnames", "chr")
  
  keep <- c("chr", "start", "end", "n_overlap_CpGs", "rid", colnames(scm))
  output <- output[, keep, with=FALSE]
  
  if(verbose) message("Region summary generating in ",stop_time())
  
  return(output)
}


#--------------------------------------------------------------------------------------------------------------------------
#' Filter matrices by coverage
#' @details Takes \code{\link{scMethrix}} object and filters CpGs based on coverage statistics
#' @inheritParams generic_scMethrix_function
#' @param cov_thr minimum coverage required to call a loci covered
#' @param min_samples Minimum number of samples that should have a loci with coverage >= \code{cov_thr}. If \code{group} is given, then this applies per group. Only need one of \code{prop_samples} or \code{min_samples}.
#' @param prop_samples Minimum proportion of samples that should have a loci with coverage >= \code{cov_thr}. If \code{group} is given, then this applies per group. Only need one of \code{prop_samples} or \code{min_samples}.
#' @param group a column name from sample annotation that defines groups. In this case, the number of min_samples will be
#' tested group-wise.
#' @importFrom methods is as new
#' @examples
#' data('scMethrix_data')
#' coverage_filter(scMethrix_data, min_sample=2)
#' @return An object of class \code{\link{scMethrix}}
#' @export
coverage_filter <- function(scm = NULL, cov_thr = 1, min_samples = NULL, prop_samples=NULL, group = NULL, n_chunks=1, n_threads=1) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (!is.null(min_samples) && !is.null(prop_samples)) {
    warning("Both min_samples and prop_samples set. Defaulting to min_samples")
    prop_samples <- NULL
  } else if (is.null(min_samples) && is.null(prop_samples)) {
    stop("Neither min_samples and prop_samples is defined")
  }
  
  if (!is.numeric(cov_thr)) stop("cov_thr is not numeric") 
  
  if (!(is.numeric(min_samples) | is.numeric(prop_samples))){
    stop("min_samples and prop_samples variables must be numeric or NULL")
  }
  
  if (!is.null(group) && !(group %in% colnames(scm@colData))){
    stop(paste("The column name ", group, " can't be found in colData. Please provid a valid group column."))
  }
  
  if (n_threads > n_chunks){
    n_chunks <- n_threads
    message("n_threads should be set to be less than or equal to n_chunks.", "\n", "n_chunks has been set 
            to be equal to n_threads = ", n_threads)
  }
  
  message("Filtering by coverage...",start_time())
  
  if (is_h5(scm)) {
    if (n_chunks == 1) {
      cov_dat = get_matrix(scm, assay = "counts")
      if (!is.null(group)) {
        stop("Groups not implemented yet")
        # row_idx <- sapply(unique(m@colData[, group]), function(c) {
        #   res <- DelayedMatrixStats::rowSums2(cov_dat[, m@colData[, 
        #                                                           group] == c] >= cov_thr, na.rm = TRUE)
        #   row_idx <- (res >= max(min_samples, ceiling(prop_samples * 
        #                                                 sum(m@colData[, group] == c))))
        # })
        # row_idx <- DelayedMatrixStats::rowAlls(row_idx)
      } else {
        res <- DelayedMatrixStats::rowSums2(cov_dat >= cov_thr, na.rm = TRUE)
        row_idx <- (res >= max(min_samples, ceiling(prop_samples * 
                                                      ncol(cov_dat))))
        # row_idx <- (res >= (if (is.null(min_samples)) min_samples else ceiling(prop_samples * 
        #  
      }
    } else {
      message("Only 1 chunk supported for now")
      # 
      # 
      # if (!is.null(group)) {
      #   row_idx <- unlist(mclapply(mc.cores = n_cores, 1:n_chunks, 
      #                              function(i) {
      #                                cov_dat = get_matrix(m[((i - 1) * ceiling(nrow(m)/n_chunks) + 1):min(i * ceiling(nrow(m)/n_chunks), nrow(m)), ], 
      #                                                     type = "C")
      #                                row_idx <- sapply(unique(m@colData[, group]), function(c) {
      #                                  res <- DelayedMatrixStats::rowSums2(cov_dat[, m@colData[,group] == c] >= cov_thr, na.rm = TRUE)
      #                                  row_idx <- (res >= max(min_samples, ceiling(prop_samples * sum(m@colData[, group] == c))))
      #                                })
      #                                row_idx <- DelayedMatrixStats::rowAlls(row_idx)
      #                              }))
      # } else {
      #   row_idx <- unlist(mclapply(mc.cores = n_cores, 1:n_chunks, 
      #                              function(i) {
      #                                cov_dat = get_matrix(m[((i - 1) * ceiling(nrow(m)/n_chunks) + 1):min(i * ceiling(nrow(m)/n_chunks), nrow(m)), ], 
      #                                                     type = "C")
      #                                res <- DelayedMatrixStats::rowSums2(cov_dat >= cov_thr, na.rm = TRUE)
      #                                row_idx <- (res >= max(min_samples, ceiling(prop_samples * ncol(cov_dat))))
      #                              }))
      # }
    }
  } else {
    cov_dat = get_matrix(scm, assay = "counts")
    if (!is.null(group)) {
      stop("Groups not implemented yet")
      # row_idx <- sapply(unique(m@colData[, group]), function(c) {
      #   res <- matrixStats::rowSums2(cov_dat[, m@colData[, group] == 
      #                                          c] >= cov_thr, na.rm = TRUE)
      #   row_idx <- (res >= max(min_samples, ceiling(prop_samples * sum(m@colData[, group] == c))))
      # })
      # row_idx <- matrixStats::rowAlls(row_idx)
    } else {
      res <- matrixStats::rowSums2(cov_dat >= cov_thr, na.rm = T)
      row_idx <- (res >= max(min_samples, ceiling(prop_samples * ncol(cov_dat))))
    }
  }
  
  gc()
  message(paste0("Retained ", format(sum(row_idx), big.mark = ","),
                 " of ", format(nrow(scm), big.mark = ","), " sites"))
  
  message("Filtered in ",stop_time())
  
  return(scm[row_idx, ])
  
}

#--------------------------------------------------------------------------------------------------------------------------
#' Saves an HDF5 \code{\link{scMethrix}} object
#' @details Takes \code{\link{scMethrix}} object and saves it in the specified directory
#' @inheritParams generic_scMethrix_function
#' @param replace Should it overwrite the pre-existing data? FALSE by default.
#' @param ... Parameters to pass to saveHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @examples
#' data('scMethrix_data')
#' dir <- paste0(tempdir(),"/h5")
#' scm <- convert_scMethrix(scMethrix_data, h5_dir=dir)
#' save_HDF5_scMethrix(scm, h5_dir = dir, replace = TRUE)
#' @return Nothing
#' @export
save_HDF5_scMethrix <- function(scm = NULL, h5_dir = NULL, replace = FALSE, ...) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is.null(dir)) {
    stop("Please provide the target directory containing ")
  }
  
  if (is_h5(scm)) {
    HDF5Array::saveHDF5SummarizedExperiment(x = scm, dir = h5_dir, replace = replace, ...)
  } else {
    stop("The object is not a methrix object or not in an HDF5 format. ")
  }
}

#--------------------------------------------------------------------------------------------------------------------------
#' Loads HDF5 \code{\link{scMethrix}} object
#' @details Takes  directory with a previously saved HDF5Array format \code{\link{scMethrix}} object and loads it
#' @param dir The directory to read in from. Default NULL
#' @param ... Parameters to pass to \code{\link{loadHDF5SummarizedExperiment}}
#' @return An object of class \code{\link{scMethrix}}
#' @examples
#' data('scMethrix_data')
#' dir <- paste0(tempdir(),"/h5")
#' scm <- convert_scMethrix(scMethrix_data, h5_dir=dir)
#' save_HDF5_scMethrix(scm, h5_dir = dir, replace = TRUE)
#' n <- load_HDF5_scMethrix(dir)
#' @export
load_HDF5_scMethrix <- function(dir = NULL, ...) {
  
  if (is.null(dir)) {
    stop("Please provide the target directory containing ")
  }
  
  message("Loading HDF5 object", start_time())
  
  scm <- HDF5Array::loadHDF5SummarizedExperiment(dir = dir, ...)
  scm <- as(scm, "scMethrix")
  
  message("Loaded in ",stop_time())
  
  return(scm)
}


#--------------------------------------------------------------------------------------------------------------------------
#' Converts HDF5 \code{\link{scMethrix}} object to an in-memory \code{\link{scMethrix}} object.
#' @details Takes an HDF%-based \code{\link{scMethrix}} object and returns with the same object with in-memory assay slots.
#' @inheritParams generic_scMethrix_function
#' @return An object of class \code{\link{scMethrix}}
#' @importFrom SummarizedExperiment assays assays<-
#' @examples
#' data('scMethrix_data')
#' dir <- paste0(tempdir(),"/h5")
#' scm <- convert_scMethrix(scMethrix_data, h5_dir=dir)
#' convert_HDF5_scMethrix(scm)
#' @export
convert_HDF5_scMethrix <- function(scm = NULL) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!is_h5(scm)) {
    stop("Input scMethrix must be in HDF5 format.")
  }
  
  message("Converting in-memory scMethrix to HDF5")#,start_time())
  
  assays(scm)[[1]] <- as.matrix(assays(scm)[[1]])
  scm@metadata$is_h5 <- FALSE
  
  # message("Converted in ",stop_time())
  
  return(scm)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Converts an in-memory \code{\link{scMethrix}} to an HDF5 \code{\link{scMethrix}}
#' @details Takes a \code{\link{scMethrix}} object and returns with the same object with delayed array assay slots
#' with HDF5 backend. Might take long time!
#' @inheritParams generic_scMethrix_function
#' @return An object of class \code{\link{scMethrix}}, HDF5 format
#' @importFrom SummarizedExperiment assays
#' @examples
#' data('scMethrix_data')
#' convert_scMethrix(scMethrix_data, h5_dir=paste0(tempdir(),"/h5"))
#' @export
convert_scMethrix <- function(scm = NULL, h5_dir = NULL, verbose = TRUE) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is_h5(scm)) {
    stop("Input scMethrix is already in HDF5 format.")
  }
  
 # if (verbose) message("Converting in-memory scMethrix to HDF5", start_time())

  scm <- create_scMethrix(assays = assays(scm), h5_dir = h5_dir,
                      rowRanges = rowRanges(scm), is_hdf5 = TRUE, genome_name = scm@metadata$genome,
                      colData = scm@colData, chrom_size = scm@metadata$chrom_sizes, 
                      desc = scm@metadata$descriptive_stats, replace = TRUE, verbose = verbose)

 # if (verbose) message("Converted in ", stop_time())

  return(scm)
}

# 
# #' Order \code{\link{scMethrix}} object by SD
# #' @details Takes \code{\link{scMethrix}} object and reorganizes the data by standard deviation
# #' @param scm \code{\link{scMethrix}} object
# #' @param zero.rm Removes zero values from equations (the default empty value for sparse matrices)
# #' @param na.rm Removes the NA values from equations
# #' @return An object of class \code{\link{scMethrix}}
# #' @export
# order_by_sd <- function (scm, zero.rm = FALSE, na.rm = FALSE) {
#   # 
#   # if (!is(m, "scMethrix")){
#   #   stop("A valid scMethrix object needs to be supplied.")
#   # }
#   # 
#   # sds <- DelayedMatrixStats::rowSds(x = get_matrix(m),na.rm=TRUE)
#   # 
#   # m$sd <- sds
#   
# }

#' Subsets an \code{\link{scMethrix}} object based on \code{regions}, \code{contigs} and/or \code{samples}.
#' @details Takes \code{\link{scMethrix}} object and filters CpGs based on region, contig and/or sample. Can 
#' either subset (\code{include}) to or filter (\code{exclude}) the specified parameters.
#' @inheritParams generic_scMethrix_function
#' @param regions genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param contigs string; array of chromosome names to subset by
#' @param samples string; array of sample names to subset by
#' @param by string to decide whether to "include" or "exclude" the given criteria from the subset
#' @importFrom IRanges subsetByOverlaps
#' @examples
#' data('scMethrix_data')
#' 
#' contigs <- c("chr1","chr3")
#' regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,100000000)) 
#' samples <- c("C1","C2")
#' 
#' #Subset to only samples bed1 and bed3, and chromosome 1
#' subset_scMethrix(scMethrix_data, samples = samples, contigs = contigs, by = "include")
#' 
#' #Subset to only region "chr1:1-5"
#' subset_scMethrix(scMethrix_data, regions = regions, by = "include")
#' 
#' #Subset to exclude samples bed1 and bed3, and chromosome 1
#' subset_scMethrix(scMethrix_data, samples = samples, contigs = contigs, by = "exclude")
#' 
#' #Subset to exclude region "chr1:1-5"
#' subset_scMethrix(scMethrix_data, regions = regions, by = "exclude")
#' @return An object of class \code{\link{scMethrix}}
#' @export
subset_scMethrix <- function(scm = NULL, regions = NULL, contigs = NULL, samples = NULL, by=c("include","exclude"), verbose=TRUE) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is.null(regions) & is.null(contigs) & is.null(samples)) {
    warning("At least 1 argument mandatory for subsetting. No subset generated")
    return(scm)
  }
  
  by <- match.arg(arg = by, choices = c("include", "exclude"), several.ok = FALSE)
  
  if (by == "exclude") {
    
    if (!is.null(regions)) {
      regions <- cast_granges(regions)
      reg <- subsetByOverlaps(rowRanges(scm), regions, invert = TRUE, type="any", maxgap=-1L, minoverlap=0L)
      scm <- subset_scMethrix(scm,regions=reduce(reg),by="include")
    }
    
    if (!is.null(contigs)) {
      c <- as.character(seqnames(scm)@values)
      scm <- subset_scMethrix(scm,contigs = c[!c %in% contigs],by="include")
    }
    
    if (!is.null(samples)) {
      s <- row.names(colData(scm))
      scm <- subset_scMethrix(scm,samples = s[!s %in% samples],by="include")
    }
    
  } else {
    
    if (verbose) message("Subsetting CpG sites...",start_time())
    
    if (!is.null(regions)) {
      regions <- cast_granges(regions)
      if (verbose) message("   Subsetting by regions")
      subset <- cast_granges(regions)
      scm <- scm[GenomicRanges::findOverlaps(rowRanges(scm), regions)@from]
      if (nrow(scm) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
    }
    
    if (!is.null(contigs)) {
      if (verbose) message("   Subsetting by contigs")
      scm <- subset(scm, subset = as.vector(seqnames(scm)) %in% contigs)
      if (nrow(scm) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
    }
    
    if (!is.null(samples)) {
      if (verbose) message("   Subsetting by samples")
      scm <- subset(scm, select = row.names(colData(scm)) %in% samples)
      if (length(scm) == 0) stop("Samples not present in the object", call. = FALSE)
    }
    
    if (verbose) message("Subset in ",stop_time())
  }
  
  return(scm)
  
}

#--------------------------------------------------------------------------------------------------------------------------
#' Estimate descriptive statistics for each sample
#' @details Calculate descriptive statistics (mean, median, SD) either by sample or \code{per_chr}
#' @inheritParams generic_scMethrix_function
#' @param per_chr Estimate stats per chromosome. Default TRUE
#' @examples
#' data('scMethrix_data')
#' 
#' #Get stats for each sample and chromosome
#' get_stats(scMethrix_data)
#' 
#' #Get stats for each sample
#' get_stats(scMethrix_data,per_chr = FALSE)
#' @return data.table of summary stats
#' @export
get_stats <- function(scm = NULL, per_chr = TRUE) {

  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  message("Getting descriptive statistics...",start_time())
  
  ends <- len <- seqnames(scm)@lengths
  for (i in 1:length(ends)) ends[i] <- sum(as.vector(len[1:i]))
  starts <- head(c(1, ends + 1), -1)
  
  if (is_h5(scm)) {
    if (per_chr) {
      stats <- lapply(1:length(starts), function(x) {
        data.table::data.table(
          Chr = levels(seqnames(scm))[x],
          Sample_Name = colnames(scm),
          mean_meth = DelayedMatrixStats::colMeans2(get_matrix(scm), rows = starts[x]:ends[x], na.rm = TRUE),
          median_meth = DelayedMatrixStats::colMedians(get_matrix(scm), rows = starts[x]:ends[x], na.rm = TRUE),
          sd_meth = DelayedMatrixStats::colSds(get_matrix(scm), rows = starts[x]:ends[x], na.rm = TRUE)
        )
      })
      
      stats <- data.table::rbindlist(l = stats, use.names = TRUE)
      
    } else {
      stats <- data.table::data.table(
        Sample_Name = colnames(scm),
        mean_meth = DelayedMatrixStats::colMeans2(get_matrix(scm), na.rm = TRUE),
        median_meth = DelayedMatrixStats::colMedians(get_matrix(scm), na.rm = TRUE),
        sd_meth = DelayedMatrixStats::colSds(get_matrix(scm), na.rm = TRUE)
      )
    }
    
  } else {
    if (per_chr) {
      stats <- lapply(1:length(starts), function(x) {
        data.table::data.table(
          Chr = levels(seqnames(scm))[x],
          Sample_Name = colnames(scm),
          mean_meth = matrixStats::colMeans2(get_matrix(scm),rows = starts[x]:ends[x],na.rm = TRUE),
          median_meth = matrixStats::colMedians(get_matrix(scm), rows = starts[x]:ends[x], na.rm = TRUE),
          sd_meth = matrixStats::colSds(get_matrix(scm),rows = starts[x]:ends[x],na.rm = TRUE)
        )
      })
      stats <- data.table::rbindlist(l = stats, use.names = TRUE)
    } else {
      stats <-
        data.table::data.table(
          Sample_Name = colnames(scm),
          mean_meth = matrixStats::colMeans2(get_matrix(scm), na.rm = TRUE),
          median_meth = matrixStats::colMedians(get_matrix(scm), na.rm = TRUE),
          sd_meth = matrixStats::colSds(get_matrix(scm), na.rm = TRUE)
        )
    }
  }
  
  gc()
  message("Finished in ", stop_time())
  
  return(stats)
}

#--------------------------------------------------------------------------------------------------------------------------

#' Extract assays from an \code{\link{scMethrix}} object
#' @details Takes \code{\link{scMethrix}} object and returns the \code{methylation} matrix
#' @inheritParams generic_scMethrix_function
#' @param add_loci Default FALSE. If TRUE adds CpG position info to the matrix and returns as a data.table
#' @param in_granges Do you want the outcome in \code{\link{GRanges}}?
#' @param order_by_sd Order output matrix by standard deviation
#' @return HDF5Matrix or matrix
#' @examples
#' data('scMethrix_data')
#' 
#' # Get methylation data
#' get_matrix(scMethrix_data)
#' 
#' # Get methylation data with loci
#' get_matrix(scMethrix_data, add_loci=TRUE)
#' 
#' # Get methylation data with loci inside a Granges object 
#' get_matrix(scMethrix_data, add_loci=TRUE, in_granges=TRUE)
#' 
#' # Get methylation data sorted by SD
#' get_matrix(scMethrix_data, order_by_sd = TRUE)
#' @export
get_matrix <- function(scm = NULL, add_loci = FALSE, in_granges=FALSE, assay = "score", order_by_sd=FALSE) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (add_loci == FALSE & in_granges == TRUE) {
    warning("Without genomic locations (add_loci= FALSE), it is not possible to convert the results to GRanges, ", 
            "the output will be a data.frame object. ")
  }
  
  assay <- match.arg(arg = assay, choices = SummarizedExperiment::assayNames(scm))
  mtx <- SummarizedExperiment::assay(x = scm, i = which(assay == SummarizedExperiment::assayNames(scm)))
  
  if (order_by_sd) {
    if (is_h5(scm)) {
      sds = DelayedMatrixStats::rowSds(mtx, na.rm = TRUE)
    } else {
      sds = matrixStats::rowSds(mtx, na.rm = TRUE)
    }
  }
  
  if (add_loci) {
    
    if (is_h5(scm)) mtx <- as.data.frame(mtx)
    
    mtx <- as.data.frame(cbind(as.data.frame(rowRanges(scm))[,1:3], mtx))
    
    if (in_granges) {
      mtx <- GenomicRanges::makeGRangesFromDataFrame(mtx, keep.extra.columns = TRUE)
      
    } else {
      data.table::setDT(x = mtx)
      colnames(mtx)[1] <- "chr" #TODO: figure out why seqnames is used instead of chr
    }
    
  }
  
  if (order_by_sd) mtx <- mtx[order(sds, decreasing = TRUE), ]
  
  return (mtx)
}

#------------------------------------------------------------------------------------------------------------

#' Remove loci that are uncovered across all samples
#' @details Takes \code{\link{scMethrix}} object and removes loci that are uncovered across all samples
#' @inheritParams generic_scMethrix_function
#' @return An object of class \code{\link{scMethrix}}
#' @examples
#' data('scMethrix_data')
#' # Remove uncovered CpGs after subsetting to a single sample
#' remove_uncovered(subset_scMethrix(scMethrix_data, samples = "df1", by="include"))
#' @export
remove_uncovered <- function(scm = NULL) {
  
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  message("Removing uncovered CpGs...", start_time())
  
  row_idx <- rowSums(!is.na(get_matrix(scm)))==0
  
  message(paste0("Removed ", format(sum(row_idx), big.mark = ","),
                 " [", round(sum(row_idx)/nrow(scm) * 100, digits = 2), "%] uncovered loci of ",
                 format(nrow(scm), big.mark = ","), " sites (",stop_time(),")"))
  
  if (!sum(row_idx) == 0) scm <- scm[!row_idx, ]
  
  return(scm)
}



#--------------------------------------------------------------------------------------------------------------------------

#' Masks high or low coverage \code{count} or \code{cell} count
#' @details Takes \code{\link{scMethrix}} object and masks sites with too high or too low coverage
#'  by putting NA for coverage and beta value. The sites will remain in the object. 
#' @inheritParams generic_scMethrix_function
#' @param low_count The minimal coverage allowed. Everything below, will get masked. Default = NULL, nothing gets masked.
#' @param high_quantile The quantile limit of coverage. Quantiles are calculated for each sample and everything that belongs to a
#' higher quantile than the defined will be masked. Default = 0.99.
#' @param n_threads Number of parallel instances. Can only be used if \code{\link{scMethrix}} is in HDF5 format. Default = 1.
#' @param type Whether to use the "counts" coverage matrix or "cells" count when masking
#' @return An object of class \code{\link{scMethrix}}
#' @importFrom SummarizedExperiment assays assays<-
#' @examples
#' data('scMethrix_data')
#' mask_scMethrix(scMethrix_data,low_count=4,type="counts")
#' @export
mask_scMethrix <- function(scm = NULL, low_count = 0, high_quantile = NULL, n_threads=1 ,type="counts") {
  
  if (!is(scm, "scMethrix")) stop("A valid scMethrix object needs to be supplied.")
  
  if (!is_h5(scm) & n_threads != 1) 
    stop("Parallel processing not supported for a non-HDF5 scMethrix object due to probable high memory usage. \nNumber of cores (n_threads) needs to be 1.")
  
  type = match.arg(arg = type, choices = c('cells', 'counts'))
  
  if (!has_cov(scm) && type == "counts") stop("No coverage matrix is present in the object. 
                                              Retry with type='cells'")
  
  if (!is.null(high_quantile)) {
    if (type == 'cells') stop("high_quantile cannot be used with 'cells' ")
    if (high_quantile >= 1 | high_quantile <= 0) stop("High quantile should be between 0 and 1. ")
  }
  
  message("Masking CpG sites in ",high_quantile," quintile and ",type, " < ",low_count, start_time())
  
  if (!is.null(low_count)) {
    
    if(!is.numeric(low_count)){
      stop("low_count must be a numeric value.")
    }
    
    n <- nrow(scm) - DelayedMatrixStats::colCounts(get_matrix(scm), value = as.integer(NA))
    
    if (type == "cells") {
      row_idx <- DelayedMatrixStats::rowCounts(get_matrix(scm), 
                                               value = as.integer(NA)) > (nrow(colData(scm))-low_count)
    } else {
      row_idx <- DelayedMatrixStats::rowSums2(get_matrix(scm,assay="counts"),na.rm=TRUE) < low_count
    }
    
    if (sum(row_idx) == 0) stop("No CpGs found with low_count")
    
    row_idx <- which(row_idx)  
    emp <- array(as.integer(NA),c(length(row_idx), nrow(colData(scm))))
    
    assays(scm)[[1]][row_idx,] <- emp
    if (type == "counts") assays(scm)[[2]][row_idx,] <- emp 
    
    n <- n-(nrow(scm) - DelayedMatrixStats::colCounts(get_matrix(scm), value = as.integer(NA)))
    
    for (i in seq_along(colnames(scm))) {
      if (n[i]==0) next
      message(paste0("   Masked ", n[i], " CpGs due to too low cells in sample ",  colnames(scm)[i], "."))
    }
  }
  
  if (!is.null(high_quantile)) {
    
    message("Masking coverage higher than ",high_quantile * 100," percentile")
    
    quantiles <- DelayedMatrixStats::colQuantiles(assays(scm)[[2]], probs = high_quantile,
                                                  na.rm = TRUE, drop = F)
    quantiles <- as.vector(quantiles)
    names(quantiles) <- rownames(scm@colData)
    
    row_idx2 <- t(t(assays(scm)[[2]]) > quantiles)
    
    if (sum(row_idx2, na.rm = TRUE) == 0) stop("No samples found within in the quantile")
    
    assays(scm)[[1]][row_idx2] <- as.integer(NA)
    assays(scm)[[2]][row_idx2] <- as.integer(NA)
    
    n <- DelayedMatrixStats::colSums2(row_idx2, na.rm = T)
    
    for (i in seq_along(colnames(scm))) {
      if (n[i]==0) next
      message(paste0("Masked ", n[i], " CpGs due to too high coverage in sample ",
                     colnames(row_idx2)[i], "."))
    }
  }
  
  message("Masked in ",stop_time())
  
  return(scm)
}