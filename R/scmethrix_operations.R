#------------------------------------------------------------------------------------------------------------
#' Adds descriptive statistics to metadata columns in an \code{\link{scMethrix}} object.
#' @details Adds the mean, median, SD, and sample count and coverage (if present) for  the \code{\link{GenomicRanges}} in an \code{\link{scMethrix}} object. This can be accessed using mcols().
#' 
#' This data will not be updated automatically for any subset, merge, bin, etc functions.
#' 
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
  
  stats <-
    data.table::data.table(
      mean_meth = DelayedMatrixStats::rowMeans2(score(scm), na.rm = TRUE),
      median_meth = DelayedMatrixStats::rowMedians(score(scm), na.rm = TRUE),
      sd_meth = DelayedMatrixStats::rowSds(score(scm), na.rm = TRUE),
      cells = ncol(scm)-DelayedMatrixStats::rowCounts(score(scm), value = NA)
    )
    
  if(has_cov(scm)) stats[,"counts" := DelayedMatrixStats::rowSums2(get_matrix(scm,assay="counts"), na.rm = TRUE)] 
  
  mcols(scm) <- stats
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Removes an assay from an \code{\link{scMethrix}} object
#' @details This will remove an assay from the scMethrix experiment object. All transformed assays may be removed, as well as the coverage assay (since it is less useful when compared to normal WGBS data), but the score assay cannot be removed. Reduced dimensionality data will be retained even if the parent assay is removed.
#' @inheritParams generic_scMethrix_function
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' remove_assay(scMethrix_data,assay="counts")
#' @export
remove_assay <- function(scm=NULL, assay=NULL) {
  
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
#' @details Merges the assay data from two \code{\link{scMethrix}} objects. Assays not shared between assays will be dropped, as well as all reduced dimensionality data.
#' 
#' Requirements for merging
#'    - If merging by rows, all CpG sites must be unique and samples must be identical
#'    - If merging by columns, all samples must be unique and CpG sites must be identical
#' 
#' Metadata will be retained in certain situations:#' 
#'    For row merges:
#'       - Ranges metadata (\code{mcols()}) will be merged, with missing columns in either assay filled with NAs
#'       - Sample metadata (\code{colData()}) will attempt to be merged. Overlapping, non-identical columns will be appended with `.1` and `.2`.
#'
#'    For row merges:
#'       - Ranges metadata (\code{mcols()}) will attempt to be merged. Overlapping, non-identical columns will be appended with `.1` and `.2`.
#'       - Sample metadata (\code{colData()}) will be merged, with missing columns in either assay filled with NAs
#' 
#'    For both merges:
#'       - Experiment metadata (\code{metadata()}) will attempt to be merged. Overlapping, non-identical columns will be appended with `.1` and `.2`.
#'  
#'    Custom experiment metadata can manually be added via \code{metadata() <-}, or to rowRanges via \code{mcols() <-}.
#' 
#' @param scm1 A \code{\link{scMethrix}} object
#' @param scm2 A \code{\link{scMethrix}} object
#' @param by Merge by columns or rows
#' @return A merged \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' merge_scMethrix(scMethrix_data[1:5],scMethrix_data[6:10],by="row")
#' merge_scMethrix(scMethrix_data[,1:2],scMethrix_data[,3:4],by="col")
#' @export
merge_scMethrix <- function(scm1 = NULL, scm2 = NULL, by = c("row", "col")) {

  if (!is(scm1, "scMethrix") || !is(scm2, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is_h5(scm1) != is_h5(scm2)) stop("Both input objects must be either in-memory or HDF5 format.")
  
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
  
  # Fix duplicate names in metadata, append if there are common elements that have different values
  slots <- c(metadata)
  
  if (by == "row") slots <- c(slots,colData,int_colData)
  if (by == "col") slots <- c(slots,mcols,elementMetadata,int_elementMetadata)

  invisible(lapply(slots, function(op) {

    if (!isTRUE(all.equal(op(scm1),op(scm2))) && length(op(scm1))>0) {

      meta1 <- intersect(names(op(scm1)), names(op(scm2)))
      meta2 <- intersect(names(op(scm2)), names(op(scm1)))
      
      if (length(c(meta1,meta2)) != 0) {
        warning("Same metadata columns are present in ",op@generic,"(). These will be appended with `.1` or `.2`")

        sapply(1:2, function(s) {
          names <- paste0("names(",op@generic,"(scm",s,"))")
          do <- paste0(names,"[",names," %in% meta",s,"] <<- paste0(",names,"[",names," %in% meta",s,"],'.",s,"')")
          eval(parse(text = eval(expression(do))))
        })
      }
    } 
  }))

  # Merge by row
  if (by == "row") {
    if (nrow(colData(scm1)) != nrow(colData(scm2)) || !all(rownames(scm1@colData) == rownames(scm2@colData))) 
      stop("You have different samples in your dataset. You need the same samples in your datasets. ")
    
    if (length(intersect(rowRanges(scm1),rowRanges(scm2))) != 0)
      stop("There are overlapping regions in your datasets. Each object must contain unique regions. ")
    
    scm <- rbind(scm1, scm2)
    scm <- sort(scm)
  }
  
  # Merge by col
  if (by == "col") {
    
    if (any(rownames(scm1@colData) %in% rownames(scm2@colData))) 
      stop("You have the same samples in your datasets. You need different samples for this merging.  ")
    
    
    if (length(intersect(rowRanges(scm1),rowRanges(scm2))) != length(rowRanges(scm1))) 
      stop("There are non-overlapping regions in your datasets. This function only takes identical regions. ")

    # Merge sample metadata. Ensure the column names match, fill with NAs if not
    colData(scm1)[setdiff(names(colData(scm2)), names( colData(scm1)))] <- NA
    colData(scm2)[setdiff(names( colData(scm1)), names(colData(scm2)))] <- NA
      
    scm <- cbind(scm1, scm2)

    # Ensure object is consistent regardless the order of scm1 and scm2
    ord <- order(colnames(scm))
    colData(scm) <- colData(scm)[ord , order(names(colData(scm))), drop=FALSE]
    
    dimnames(scm)[[2]] <- sort(colnames(scm))
    
    for (name in SummarizedExperiment::assayNames(scm)) {
      assay(scm,name,withDimnames = F) <- assay(scm,name,withDimnames = F)[ , ord]
    }
  }
  
  #Realize if HDF5
  if (is_h5(scm)) {
    for (name in SummarizedExperiment::assayNames(scm)) {
       assay(scm,name) <- as(assay(scm,name), "HDF5Array")
    }
  }
  
  # Remove duplicate experiment metadata
  invisible(lapply(c(metadata,int_metadata), function(op) {
    eval(parse(text = eval(expression(paste0(op@generic,"(scm) <<- ",op@generic,"(scm)[unique(names(",op@generic,"(scm)))]")))))
  }))

  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Extracts and summarizes methylation or coverage info by regions of interest
#' @details Takes \code{\link{scMethrix}} object and summarizes regions
#' @inheritParams generic_scMethrix_function
#' @param regions genomic regions to be summarized. Could be a data.table with 3 columns (chr, start, end) or a \code{\link{GenomicRanges}} object
#' @param type string; the assay to be summarized. Default 'score'
#' @param how closure; mathematical function by which regions should be summarized. Can be one of the following: mean, sum, max, min. Default 'mean'
#' @param overlap_type string; defines the type of the overlap of the CpG sites with the target region. Default value is `within`. For detailed description, see the \code{findOverlaps} function of the \code{\link{IRanges}} package.
#' #param elementMetadata.col columns in \code{\link{scMethrix}}@elementMetadata which needs to be summarised. Default = NULL.
#' @param n_chunks Number of chunks to split the \code{\link{scMethrix}} object in case it is very large. Default = 1.
#' @param n_threads Number of parallel instances. \code{n_cores} should be less than or equal to \code{n_chunks}. If \code{n_chunks} is not specified, then \code{n_chunks} is initialized to be equal to \code{n_cores}. Default = 1.
#' @param group a column name from sample annotation that defines groups. In this case, the number of min_samples will be tested group-wise.
#' @importFrom methods setClass
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
  
  if (!is.null(regions)) {
    regions = cast_granges(regions)
  } else { # If no region is specifed, use entire chromosomes
    regions = range(rowRanges(scm))
  }
  
  regions$rid <- paste0("rid_", 1:length(regions))
  
  overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, regions, type = overlap_type)) #GenomicRanges::findOverlaps(rowRanges(m), regions)@from
  
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
    
    #stop("Chunking not enabled")

    if(nrow(overlap_indices) < n_chunks){
      n_chunks <- nrow(overlap_indices)
      warning("Fewer overlaps indicies than n_chunks. Defaulting to n_chunks = ",n_chunks)
    }

    if (n_chunks < n_threads) {
      n_threads <- n_chunks
      warning("n_threads < n_chunks. Defaulting to n_threads = ",n_threads)
    }

    cl <- parallel::makeCluster(n_threads)
    doParallel::registerDoParallel(cl)

    parallel::clusterEvalQ(cl, c(library(data.table), library(SingleCellExperiment), sink(paste0("D:/Git/scMethrix/", Sys.getpid(), ".txt"))))
    parallel::clusterEvalQ(cl, expr={
      scMethrix <- setClass(Class = "scMethrix", contains = "SingleCellExperiment")
    })
    parallel::clusterExport(cl,list('scm','scMethrix','type','is_h5','get_matrix','start_time','split_time','stop_time'))

    chunk_overlaps <- split(overlap_indices$xid, ceiling(seq_along(overlap_indices$xid) /
                                                        ceiling(length(overlap_indices$xid)/n_chunks)))

    data <- parallel::parLapply(cl,chunk_overlaps,fun=function(i) {
      get_matrix(scm[i,], assay = type, add_loci = TRUE) # TODO: object of type 'S4' is not subsettable
    })

    parallel::stopCluster(cl)

    data <- rbindlist(data)
    
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
#' Extract assays from an \code{\link{scMethrix}} object
#' @details Takes \code{\link{scMethrix}} object and returns the \code{methylation} matrix. This will return in the format used by the object (matrix or HDF5matrix).
#' @inheritParams generic_scMethrix_function
#' @param add_loci Default FALSE. If TRUE adds CpG position info to the matrix and returns as a data.table
#' @param in_granges Do you want the outcome in \code{\link{GRanges}}?
#' @param order_by_sd Order output matrix by standard deviation
#' @return HDF5Matrix or matrix
#' @import SummarizedExperiment
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
#--------------------------------------------------------------------------------------------------------------------------
get_matrix <- function(scm = NULL, assay = "score", add_loci = FALSE, in_granges=FALSE, order_by_sd=FALSE) {
  
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
save_HDF5_scMethrix <- function(scm = NULL, h5_dir = NULL, replace = FALSE, verbose = TRUE, ...) {

  if (is(scm, "scMethrix")) {
    if (!is_h5(scm)) stop("A valid scMethrix HDF5 object needs to be supplied.", call. = FALSE)
  } else if (!is(scm, "SingleCellExperiment")) {
    stop("A valid SingleCellExperiment-derived object needs to be supplied.", call. = FALSE)
  }

  if (is.null(h5_dir)) {
    stop("Please provide the target directory containing ")
  }

  #if (is.null(h5_dir)) h5_dir = paste0(tempdir(),"/h5")
  
  if (dir.exists(h5_dir)) {
    files <- list.files (h5_dir,full.names = TRUE)
    
    if (length(files) != 0 && replace == FALSE) {
      message("Files are present in the target directory, including: ")
      writeLines(paste("   ",head(files)))
      choice <- menu(c("Yes", "No"), title="Are you sure you want to delete this directory and save the HDF5 experiment?")
  
      if (choice == 2 || choice == 0) {
        message("Saving aborted. The target directory has not been affected.")
        return(invisible(NULL))
      } else {
        #unlink(h5_dir, recursive=TRUE)
        replace = TRUE
      }
    }
  }
  
  if (verbose) message("Saving HDF5 experiment to disk...",start_time())

  HDF5Array::saveHDF5SummarizedExperiment(x = scm, dir = h5_dir, replace = replace, chunkdim = c(length(rowRanges(scm)),1), ...)

  if (verbose) message("Experiment saved in ",stop_time())
}

#' Loads HDF5 \code{\link{scMethrix}} object
#' @details Takes  directory with a previously saved HDF5Array format \code{\link{scMethrix}} object and loads it
#' @inheritParams generic_scMethrix_function
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
load_HDF5_scMethrix <- function(dir = NULL, verbose = TRUE, ...) {
  
  if (is.null(dir)) {
    stop("Please provide the target directory containing ")
  }
  
  if (verbose) message("Loading HDF5 object", start_time())
  
  scm <- HDF5Array::loadHDF5SummarizedExperiment(dir = dir, ...)
  scm <- as(scm, "scMethrix")
  
  if (verbose) message("Loaded in ",stop_time())
  
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
#' scm <- convert_scMethrix(scMethrix_data)
#' convert_HDF5_scMethrix(scm)
#' @export
convert_HDF5_scMethrix <- function(scm = NULL, verbose = TRUE) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (!is_h5(scm)) {
    stop("Input scMethrix must be in HDF5 format.")
  }
  
  if (verbose) message("Converting HDF5 scMethrix to in-memory", start_time())
  
  assays(scm)[[1]] <- as.matrix(assays(scm)[[1]])
  scm@metadata$is_h5 <- FALSE
  
  if (verbose) message("Converted in ", stop_time())
  
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
  
  if (verbose) message("Converting in-memory scMethrix to HDF5", start_time())

  scm <- create_scMethrix(assays = assays(scm), h5_dir = h5_dir,
                      rowRanges = rowRanges(scm), is_hdf5 = TRUE, genome_name = scm@metadata$genome,
                      colData = scm@colData, chrom_size = scm@metadata$chrom_sizes, 
                      desc = scm@metadata$descriptive_stats, replace = TRUE, verbose = verbose)

  if (verbose) message("Converted in ", stop_time())

  return(scm)
}

#' Subsets an \code{\link{scMethrix}} object based on \code{regions}, \code{contigs} and/or \code{samples}.
#' @details Takes \code{\link{scMethrix}} object and filters CpGs based on region, contig and/or sample. Can 
#' either subset (\code{include}) to or filter (\code{exclude}) the specified parameters.
#' @inheritParams generic_scMethrix_function
#' @param regions genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
#' @param contigs string; array of chromosome names to subset by
#' @param samples string; array of sample names to subset by
#' @param overlap_type string; defines the type of the overlap of the CpG sites with the target region. Default value is `within`. For detailed description, see the \code{findOverlaps} function of the \code{\link{IRanges}} package.
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
subset_scMethrix <- function(scm = NULL, regions = NULL, contigs = NULL, samples = NULL, by=c("include","exclude"), overlap_type="within",verbose=TRUE) {
  
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
      reg <- subsetByOverlaps(rowRanges(scm), regions, invert = TRUE, type=overlap_type, maxgap=-1L, minoverlap=0L)
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
#' @param per_chr boolean; Estimate stats per chromosome. Default TRUE
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
get_stats <- function(scm = NULL, assay="score",per_chr = TRUE) {

  x <- NULL
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  message("Getting descriptive statistics...",start_time())

  ends <- len <- seqnames(scm)@lengths
  for (i in 1:length(ends)) ends[i] <- sum(as.vector(len[1:i]))
  starts <- head(c(1, ends + 1), -1)
  
    if (per_chr) {
      stats <- lapply(1:length(starts), function(x) {
        s <- data.table::data.table(
          Chr = levels(seqnames(scm))[x],
          Sample_Name = colnames(scm),
          mean_meth = DelayedMatrixStats::colMeans2(get_matrix(scm,assay=assay), rows = starts[x]:ends[x], na.rm = TRUE),
          median_meth = DelayedMatrixStats::colMedians(get_matrix(scm,assay=assay), rows = starts[x]:ends[x], na.rm = TRUE),
          sd_meth = DelayedMatrixStats::colSds(get_matrix(scm,assay=assay), rows = starts[x]:ends[x], na.rm = TRUE)
        )
      })
      
      stats <- data.table::rbindlist(l = stats, use.names = TRUE)
      
    } else {
      stats <- data.table::data.table(
        Sample_Name = colnames(scm),
        mean_meth = DelayedMatrixStats::colMeans2(get_matrix(scm,assay=assay), na.rm = TRUE),
        median_meth = DelayedMatrixStats::colMedians(get_matrix(scm,assay=assay), na.rm = TRUE),
        sd_meth = DelayedMatrixStats::colSds(get_matrix(scm,assay=assay), na.rm = TRUE)
      )
    }

  gc()
  message("Finished in ", stop_time())
  
  return(stats)
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
remove_uncovered <- function(scm = NULL, verbose = TRUE) {
  
  if (!is(scm, "scMethrix")){
    stop("A valid scMethrix object needs to be supplied.")
  }
  
  if (verbose) message("Removing uncovered CpGs...", start_time())
  
  row_idx <- DelayedMatrixStats::rowSums2(!is.na(get_matrix(scm)))==0
  
  if (verbose) message(paste0("Removed ", format(sum(row_idx), big.mark = ","),
                 " [", round(sum(row_idx)/nrow(scm) * 100, digits = 2), "%] uncovered loci of ",
                 format(nrow(scm), big.mark = ","), " sites (",stop_time(),")"))
  
  if (!sum(row_idx) == 0) scm <- scm[!row_idx, ]
  
  return(scm)
}



#--------------------------------------------------------------------------------------------------------------------------
#' Masks CpGs by coverage
#' @details Takes \code{\link{scMethrix}} object and masks sites with low overall or high average coverage by putting NA for assay values. The sites will remain in the object and all assays will be affected.
#'  
#'  \code{low_threshold} is used to mask sites with low overall coverage.
#'  \code{avg_threshold} is used to mask sites with high aberrant counts. For single cell data, this is typically CpG sites with an average count > 2, as there are only two strands in a cell to sequence.
#'  
#' @inheritParams generic_scMethrix_function
#' @param low_threshold numeric; The minimal coverage allowed. Everything below will get masked. If NULL, this will be ignored.
#' @param avg_threshold numeric; The max average coverage. Typical value is 2, as there should only be 2 reads per cell. If NULL, this will be ignored. 
#' @param n_threads integer; Number of parallel instances. Can only be used if \code{\link{scMethrix}} is in HDF5 format. Default = 1
#' @return An object of class \code{\link{scMethrix}}
#' @importFrom SummarizedExperiment assays assays<-
#' @examples
#' data('scMethrix_data')
#' mask_by_coverage(scMethrix_data,low_threshold=2, avg_threshold=2)
#' @export
mask_by_coverage <- function(scm = NULL, assay = "score", low_threshold = NULL, avg_threshold = NULL, n_threads=1 , verbose = TRUE) {
  if (!is(scm, "scMethrix")) stop("A valid scMethrix object needs to be supplied.")
  
  if (!is_h5(scm) && n_threads != 1) 
    stop("Parallel processing not supported for a non-HDF5 scMethrix object due to probable high memory usage. \nNumber of cores (n_threads) needs to be 1.")
  
  if (!is.null(low_threshold) && (!is.numeric(low_threshold) || low_threshold < 0)) stop("low_threshold must be greater than 0")
  if (!is.null(avg_threshold) && (!is.numeric(avg_threshold) || avg_threshold < 0)) stop("avg_threshold must be greater than 0")
  
  if (!has_cov(scm)) stop("Cannot mask as no coverage matrix is present in the object.")
  
  # if(!is.numeric(low_threshold) || !is.numeric(low_threshold)){
  #   stop("Thresholds must be a numeric value.")
  # }
 
  if (verbose) message("Masking by coverage...",start_time())
  
  avg_row_idx <- low_row_idx <- NULL
  
  if (!is.null(low_threshold)) {
    
    if (verbose) message("Finding low coverage CpG sites...")
    
    low_row_idx <- which(DelayedMatrixStats::rowSums2(get_matrix(scm,assay="counts"),na.rm=TRUE) < low_threshold)

    if (verbose) message("   Found ",length(low_row_idx)," CpGs with coverage < ", low_threshold)
    
  }
  
  if (!is.null(avg_threshold)) {
    
    if (verbose) message("Finding high average count CpG sites...")
    
    avg_row_idx <- which(DelayedMatrixStats::rowMeans2(counts(scm), na.rm = TRUE) > avg_threshold)
    
    if (verbose) message("   Found ",length(avg_row_idx)," CpGs with average coverage > ",avg_threshold)
    
  }
  
  row_idx <- unique(c(low_row_idx,avg_row_idx))
  
  scm <- mask_helper(scm, row_idx, verbose = verbose)

  return(scm)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Masks CpGs by cell count
#' @details Takes \code{\link{scMethrix}} object and masks sites with too high or too low coverage
#'  by putting NA for assay values. The sites will remain in the object and all assays will be affected.
#'  
#'  \code{low_threshold} is used to mask sites with low overall cell counts. A site represented by a single sample is typically not useful.
#'  \code{prop_threshold} is used to mask sites with a low proportional count  
#'  
#' @inheritParams generic_scMethrix_function
#' @param low_threshold numeric; The minimal cell count allowed. Everything below will get masked.
#' @param prop_threshold numeric; The minimal proportion of covered cells.
#' @param n_threads integer; Number of parallel instances. Can only be used if \code{\link{scMethrix}} is in HDF5 format. Default = 1.
#' @return An object of class \code{\link{scMethrix}}
#' @importFrom SummarizedExperiment assays assays<-
#' @examples
#' data('scMethrix_data')
#' mask_by_sample(scMethrix_data,low_threshold=2)
#' mask_by_sample(scMethrix_data,prop_threshold=0.5) 
#' @export
mask_by_sample <- function(scm = NULL, assay = "score", low_threshold = NULL, prop_threshold = NULL, n_threads=1 , verbose = TRUE) {
  
  if (!is(scm, "scMethrix")) stop("A valid scMethrix object needs to be supplied.")
  
  if (!is_h5(scm) & n_threads != 1) 
    stop("Parallel processing not supported for a non-HDF5 scMethrix object due to probable high memory usage. \nNumber of cores (n_threads) needs to be 1.")
  
  if (!is.null(low_threshold) && !is.null(prop_threshold))
    stop("low_threshold and prop_threshold are mutually exclusive. Once must be set to NULL")
  
  # if(!is.numeric(low_threshold) || !is.numeric(prop_threshold)){
  #   stop("Thresholds must be a numeric value.")
  # }
  
  if (!is.null(low_threshold) && (!is.numeric(low_threshold) || low_threshold < 0)) stop("low_threshold must be >= 0")
  if (!is.null(prop_threshold) && (!is.numeric(prop_threshold) || prop_threshold > 1 || prop_threshold < 0)) stop("prop_threshold must be between 0 and 1")
  
  if (verbose) message("Masking by sample count...",start_time())
  
  row_idx <- NULL
  
  if (!is.null(low_threshold)) {
    
    if (verbose) message("Finding low coverage CpG sites...")
    
    row_idx <- which(DelayedMatrixStats::rowCounts(get_matrix(scm), 
                                                   value = as.integer(NA)) > (ncol(scm)-low_threshold))
    
    if (verbose) message("   Found ",length(row_idx)," CpGs with count < ",low_threshold)
    
  } else if (!is.null(prop_threshold)) {
    
    if (verbose) message("Finding high average count CpG sites...")
    
    row_idx <- which(DelayedMatrixStats::rowCounts(get_matrix(scm), 
                                                   value = as.integer(NA)) > ncol(scm)*(1-prop_threshold))
    
    if (verbose) message("   Found ",length(row_idx)," CpGs with missing cell proportion < ",prop_threshold)
    
  }
  
  scm <- mask_helper(scm, row_idx, verbose = verbose)
  
  return(scm)
}

#' Masks non-variable CpGs
#' @details Takes \code{\link{scMethrix}} object and masks CpGs with low variability by putting NA for assay values. The sites will remain in the object and all assays will be affected.
#' @inheritParams generic_scMethrix_function
#' @param low_threshold numeric; The variability threshold. Masking is done for all CpGs less than or equal to this value. A CpG that is methylated or unmethylated in all samples will have a variability of 0, whereas a completely variable CpG will have a value of 1. Default = 0.05.
#' @param n_threads integer; Number of parallel instances. Can only be used if \code{\link{scMethrix}} is in HDF5 format. Default = 1.
#' @return An object of class \code{\link{scMethrix}}
#' @importFrom SummarizedExperiment assays assays<-
#' @examples
#' data('scMethrix_data')
#' mask_by_variance(scMethrix_data,low_threshold=0.05)
#' @export
mask_by_variance <- function(scm = NULL, assay = "score", low_threshold = 0.05, n_threads = 1, verbose = TRUE) {
  
  if (!is(scm, "scMethrix")) stop("A valid scMethrix object needs to be supplied.")
  
  if (!is_h5(scm) & n_threads != 1) 
    stop("Parallel processing not supported for a non-HDF5 scMethrix object due to probable high memory usage. \nNumber of cores (n_threads) needs to be 1.")
  
  if (!is.numeric(low_threshold) || low_threshold > 1 || low_threshold < 0) stop("low_threshold must be between 0 and 1")
  
  if (verbose) message("Masking non-variable CpG sites...",start_time())

  n <- DelayedMatrixStats::colCounts(get_matrix(scm), value = as.integer(NA))
    
  vals <- DelayedMatrixStats::rowSds(get_matrix(scm,assay=assay), na.rm = TRUE) 
  vals[which(is.na(vals))] <- 0
  vals <- vals/DelayedMatrixStats::rowMeans2(get_matrix(scm,assay=assay), na.rm = TRUE)
  vals[which(is.nan(vals))] = 0 #For case of all non-methylated CpGs (hence 0/0 = NaN)
  
  row_idx <- which(vals <= low_threshold)
    
  scm <- mask_helper(scm, row_idx, verbose = verbose)

  return(scm)
  
}

#' Helper function for masking. All rows in row_idx will be set to NA.
#' @details Used with mask_by_sample, mask_by_coverage, and mask_by_variance. It iterates through all assays in the inputted \code{\link{scMethrix}} object and replaces all rows in row_idx with NA.  
#' @inheritParams generic_scMethrix_function
#' @param row_idx numeric; A vector of row indexes for which to replace with NA
#' @return An object of class \code{\link{scMethrix}}
#' @importFrom SummarizedExperiment assays assays<-
#' @export
mask_helper <- function (scm, row_idx, verbose = TRUE) {
  
  if (length(row_idx) == 0) {
    if (verbose) message("   No CpGs were masked (",stop_time(),"elapsed)")
  } else if (length(row_idx) == nrow(scm)) {
    stop("All CpGs were masked. Re-check threshold values.")
  } else {
    for (i in 1:length(assays(scm))) {
      assays(scm)[[i]][row_idx,] <- as.integer(NA)
    } 
    
    if (verbose) message("Masked ",length(row_idx), 
                         " [", round(length(row_idx)/nrow(scm) * 100, digits = 2), "%] CpG sites in ",stop_time())
  }
  
  return(scm)

}

