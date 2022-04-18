
#---- add_assay -------------------------------------------------------------------------------------------
#' Adds an assay from an [`scMethrix`] object
#' @description Fulfills the same function as `assay(scm, assay) <- matrix`, but with additional checks.
#' @inheritParams generic_scMethrix_function
#' @param matrix `matrix`; the input matrix
#' @return An [`scMethrix`] object
#' @examples
#' data('scMethrix_data')
#' @export
add_assay <- function(scm=NULL, new_assay ="new_assay", matrix=NULL) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  
  if (!.validateAssay(scm,new_assay,is.absent=T))
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)
  
  #---- Function code ------------------------------------------------------
  assay(scm, assay) <- matrix
  
  validObject(scm)
  return(scm)
}

#---- remove_assay -------------------------------------------------------------------------------------------
#' Removes an assay from an [`scMethrix`] object
#' @description This will remove an assay from the experiment object. All transformed assays may be removed, as well as the coverage assay (since it is less useful when compared to normal WGBS data), but the score assay cannot be removed. Reduced dimensionality data will be retained even if the parent assay is removed.
#' @inheritParams generic_scMethrix_function
#' @return An [`scMethrix`] object
#' @examples
#' data('scMethrix_data')
#' remove_assay(scMethrix_data,assay="counts")
#' @export
remove_assay <- function(scm=NULL, assay=NULL) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  if (assay == "score") stop("Score assay cannot be removed.", call. = FALSE)
  
  #---- Function code ------------------------------------------------------
  assays(scm) <- SummarizedExperiment::assays(scm)[-which(SummarizedExperiment::assayNames(scm) == assay)]
  
  validObject(scm)
  return(scm)
}

#---- merge_scMethrix --------------------------------------------------------------------------------------------------
#' Merges two [`scMethrix-class`] objects by `row` or `col`
#' @details Merges the assay data from two [`scMethrix-class`] objects. Assays not shared between assays will be dropped, as well as all reduced dimensionality data.
#' 
#' Requirements for merging
#'    - If merging by rows, all CpG sites must be unique and samples must be identical
#'    - If merging by columns, all samples must be unique and CpG sites must be identical
#' 
#' Metadata will be retained in certain situations:
#'    For row merges:
#'       - Ranges metadata (`mcols()`) will be merged, with missing columns in either assay filled with NAs
#'       - Sample metadata (`colData()`) will attempt to be merged. Overlapping, non-identical columns will be appended with `.1` and `.2`.
#'
#'    For row merges:
#'       - Ranges metadata (`mcols()`) will attempt to be merged. Overlapping, non-identical columns will be appended with `.1` and `.2`.
#'       - Sample metadata (`colData()`) will be merged, with missing columns in either assay filled with NAs
#' 
#'    For both merges:
#'       - Experiment metadata (`metadata()`) will attempt to be merged. Overlapping, non-identical elements will be appended with `.1` and `.2`.
#'  
#'    Custom experiment metadata can manually be added via `metadata() <-`, or to rowRanges via `mcols() <-`.
#' @inheritParams generic_scMethrix_function
#' @param scm1 [`scMethrix-class`]; A single cell methylation experiment
#' @param scm2 [`scMethrix-class`]; A single cell methylation experiment
#' @param by `string`; Merge by 'columns' or 'rows'
#' @return A merged [`scMethrix-class`] object
#' @examples
#' data('scMethrix_data')
#' merge_scMethrix(scMethrix_data[1:5],scMethrix_data[6:10],by="row")
#' merge_scMethrix(scMethrix_data[,1:2],scMethrix_data[,3:4],by="col")
#' @export
merge_scMethrix <- function(scm1 = NULL, scm2 = NULL, h5_dir = NULL, by = c("row", "column"), verbose = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm1)
  .validateExp(scm2)
  .validateType(verbose,"boolean")
  .validateType(h5_dir,c("string","null"))
  by <- .validateArg(by,merge_scMethrix)
  
  if (is_h5(scm1) != is_h5(scm2)) stop("Both input objects must be either in-memory or HDF5 format.", call. = FALSE)
  #TODO: Not sure if above check is needed
  
  #---- Function code ------------------------------------------------------
  names1 = SummarizedExperiment::assayNames(scm1)
  names2 = SummarizedExperiment::assayNames(scm2)
  
  if (verbose) message("Merging experiment metadata")
  
  if (!all((sort(names1)==sort(names2)))) {
    warning("Assay list not identical. All non-identical assays will be dropped from merged object.")
    a1 <- intersect(names1, names2)
    a2 <- intersect(names2, names1)
    SummarizedExperiment::assays(scm1) <- SummarizedExperiment::assays(scm1)[a1]
    SummarizedExperiment::assays(scm2) <- SummarizedExperiment::assays(scm2)[a2]
  } 
  
  if (verbose) message("Merging assays...")
  
  # Merge the rest of the experiment metadata
  if (by == "row") {
    slots <- c(S4Vectors::metadata,colData)
  }
  else {
    slots <- c(S4Vectors::metadata,rowData)
  }
  
  # slots <- list(S4Vectors::metadata)
  
  for (i in 1:length(slots)) {
    
    op <- slots[[i]]
    
    op1 <- op(scm1)
    op2 <- op(scm2)
    n1 <- names(op1)
    n2 <- names(op2)
    meta <- c(op1[setdiff(n1, n2)],op2[setdiff(n2, n1)])
    not_shown = T
    
    for (n in intersect(n1,n2)) {
      if (identical(op1[n],op2[n]) || is.null(unlist(op2[n]))) {
        meta <- append(meta,op1[n])
      } else if (is.null(unlist(op1[n]))) {
        meta <- append(meta,op2[n])
      } else {
        if(not_shown) {warning("Same metadata columns are present in ",op@generic,
                               "(). These will be appended with `.1` or `.2`")}
        meta[[paste0(n,".1")]] <- unname(unlist(op1[n]))
        meta[[paste0(n,".2")]] <- unname(unlist(op2[n]))
        not_shown = F
      }
    }
    
    eval(parse(text=eval(expression(paste0(op@generic,"(scm1) <- meta")))))
    blank <- meta[-(1:length(names(meta)))]
    eval(parse(text=eval(expression(paste0(op@generic,"(scm2) <- blank")))))
    
  }
  # invisible(lapply(slots, function(op) {
  #   
  #  
  # }))
  
  # Merge by row
  if (by == "row") {
    if (nrow(SummarizedExperiment::colData(scm1)) != nrow(SummarizedExperiment::colData(scm2)) || 
        !all(rownames(scm1@colData) == rownames(scm2@colData))) 
      stop("You have different samples in your dataset. You need the same samples in your datasets.", call. = FALSE)
    
    if (length(intersect(ranges(scm1),ranges(scm2))) != 0)
      stop("There are overlapping regions in your datasets. Each object must contain unique regions.", call. = FALSE)
    
    scm <- rbind(scm1, scm2)
    scm <- sort(scm)
  } else {
    
    # Merge by col
    if (any(rownames(scm1@colData) %in% rownames(scm2@colData))) 
      stop("You have the same samples in your datasets. You need different samples for this merging.", call. = FALSE)
    
    
    #if (length(intersect(SummarizedExperiment::rowRanges(scm1),SummarizedExperiment::rowRanges(scm2))) != 
    #    length(SummarizedExperiment::rowRanges(scm1))) 
    
    if (!identical(ranges(scm1),ranges(scm2)))  {
      stop("There are non-overlapping regions in your datasets. This function only takes identical regions.", call. = FALSE)
    }
    
    # Merge sample metadata. Ensure the column names match, fill with NAs if not
    colData(scm1)[setdiff(names(SummarizedExperiment::colData(scm2)), names(SummarizedExperiment::colData(scm1)))] <- NA
    colData(scm2)[setdiff(names(SummarizedExperiment::colData(scm1)), names(SummarizedExperiment::colData(scm2)))] <- NA
    
    scm <- cbind(scm1, scm2)
    
    # Ensure object is consistent regardless the order of scm1 and scm2
    ord <- order(colnames(scm))
    SummarizedExperiment::colData(scm) <- SummarizedExperiment::colData(scm)[ord ,
                                                                             order(names(SummarizedExperiment::colData(scm))), drop=FALSE]
    
    dimnames(scm)[[2]] <- sort(colnames(scm))
    
    for (name in SummarizedExperiment::assayNames(scm)) {
      SummarizedExperiment::assay(scm,name,withDimnames = F) <- 
        SummarizedExperiment::assay(scm,name,withDimnames = F)[ , ord]
    }
  }
  
  #Realize if HDF5
  # if (is_h5(scm)) {
  #   if (verbose) message("Realizing assays to HDF5...")
  #   
  #   for (name in SummarizedExperiment::assayNames(scm)) {
  #     SummarizedExperiment::assay(scm,name) <- as(SummarizedExperiment::assay(scm,name), "HDF5Array")
  #   }
  # }
  
  validObject(scm)
  return(scm)
}


#' Does merge
#'
#' @param scm1 first scm
#' @param scm2 second scm
#' @param h5_dir dir
#' @param verbose boolean
#' @param by_row_name boolean; uses row name instead of start position
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{# TODO: add example }
merge_scMethrix2 <- function(scm1 = NULL, scm2 = NULL, h5_dir = NULL, by_row_name = FALSE ,verbose = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm1)
  .validateExp(scm2)
  .validateType(verbose,"boolean")
  
  if (is_h5(scm1) != is_h5(scm2)) stop("Both input objects must be either in-memory or HDF5 format.", call. = FALSE)
  #TODO: Not sure if above check is needed
  
  if (any(sampleNames(scm1) %in% sampleNames(scm2))) 
    stop("Experiments must contain unique set of sample names.", call. = FALSE)
  
  #---- Function code ------------------------------------------------------
  if (verbose) message("Merging assays...")
  
  names1 = SummarizedExperiment::assayNames(scm1)
  names2 = SummarizedExperiment::assayNames(scm2)
  
  # Check for dissimilar assay list
  if (!all((sort(names1)==sort(names2)))) {
    warning("Assay list not identical. All non-identical assays will be dropped from merged object.")
    a1 <- intersect(names1, names2)
    a2 <- intersect(names2, names1)
    SummarizedExperiment::assays(scm1) <- SummarizedExperiment::assays(scm1)[a1]
    SummarizedExperiment::assays(scm2) <- SummarizedExperiment::assays(scm2)[a2]
  }
  
  if (by_row_name) {
    
    #add.gr.scm1 <- setdiff(names(rowRanges(scm2)),names(rowRanges(scm1)))
    #add.gr.scm2 <- setdiff(names(rowRanges(scm1)),names(rowRanges(scm2)))
    add.gr.scm1 <- rowRanges(scm2)[which(!(names(rowRanges(scm2)) %in% names(rowRanges(scm1))))]
    add.gr.scm2<- rowRanges(scm1)[which(!(names(rowRanges(scm1)) %in% names(rowRanges(scm2))))]
    
  } else {
    
    # Add in missing rowRanges
    gr <- c(rowRanges(scm1),rowRanges(scm2))
    ##add.gr.scm1 <- gr[-S4Vectors::queryHits(findOverlaps(gr, rowRanges(scm1), type="any")),]
    ##add.gr.scm2 <- gr[-S4Vectors::queryHits(findOverlaps(gr, rowRanges(scm2), type="any")),]
    add.gr.scm1 <- rowRanges(scm2)[-S4Vectors::queryHits(findOverlaps(rowRanges(scm2), rowRanges(scm1), type="any")),]
    add.gr.scm2 <- rowRanges(scm1)[-S4Vectors::queryHits(findOverlaps(rowRanges(scm1), rowRanges(scm2), type="any")),]
  }
  
  if (length(add.gr.scm1) > 0) {
    mtx <- matrix(nrow=length(add.gr.scm1), ncol = ncol(scm1))
    assays = sapply(assayNames(scm1),function(x) mtx,simplify = FALSE,USE.NAMES = TRUE)
    scm.temp <- scMethrix(assays = assays,rowRanges = add.gr.scm1, colData <- colData(scm1))
    scm1 <- rbind(scm1,scm.temp)
  }
  
  if (length(add.gr.scm2) > 0) {
    mtx <- matrix(nrow=length(add.gr.scm2), ncol = ncol(scm2))
    assays = sapply(assayNames(scm2),function(x) mtx,simplify = FALSE,USE.NAMES = TRUE)
    scm.temp <- scMethrix(assays = assays,rowRanges = add.gr.scm2, colData <- colData(scm2))
    scm2 <- rbind(scm2,scm.temp)
  }
  
  GenomeInfoDb::seqlevels(scm2) <- GenomeInfoDb::seqlevels(scm1) #TODO: Not sure why this is necessary
  
  scm1 <- sort(scm1)
  scm2 <- sort(scm2)

  stopifnot(length(rowRanges(scm1)) == length(rowRanges(scm2)))
  
  # Combine all metadata except colData
  slots <- c(S4Vectors::metadata, SummarizedExperiment::rowData)
  
  for (i in 1:length(slots)) {
    
    op <- slots[[i]]
    
    op1 <- op(scm1)
    op2 <- op(scm2)
    n1 <- names(op1)
    n2 <- names(op2)
    meta <- c(op1[setdiff(n1, n2)],op2[setdiff(n2, n1)])
    show_warning = F
    
    for (n in intersect(n1,n2)) {
      if (identical(op1[n],op2[n]) || is.null(unlist(op2[n]))) {
        meta <- append(meta,op1[n])
      } else if (is.null(unlist(op1[n]))) {
        meta <- append(meta,op2[n])
      } else {
        show_warning = TRUE
        meta[[paste0(n,".1")]] <- unname(unlist(op1[n]))
        meta[[paste0(n,".2")]] <- unname(unlist(op2[n]))
      }
    }
    
    if(show_warning) {warning("Same metadata columns are present in ",op@generic,
                              "(). These will be appended with `.1` or `.2`")}    
    
    eval(parse(text=eval(expression(paste0(op@generic,"(scm1) <- meta")))))
    blank <- meta[-(1:length(names(meta)))]
    eval(parse(text=eval(expression(paste0(op@generic,"(scm2) <- blank")))))
  }
  
  S4Vectors::metadata(scm2)$is_h5 <- is_h5(scm1)
  
  # Combine rest of metadata
  colData(scm1)[setdiff(names(colData(scm2)), names(colData(scm1)))] <- NA
  colData(scm2)[setdiff(names(colData(scm1)), names(colData(scm2)))] <- NA
  
  scm <- cbind(scm1,scm2)
  
  return(scm)
}

#---- get_matrix -------------------------------------------------------------------------------------------------------
#' Extract assays from an [`scMethrix-class`] object
#' @description Takes an [`scMethrix-class`] object and returns an assay in a specified matrix in the format used by the object (`matrix` or `HDF5matrix`). 
#' @inheritParams generic_scMethrix_function
#' @param add_loci `boolean`; Adds genomic loci to the output. Default = `FALSE`. If `TRUE`, it adds CpG position info to the matrix and returns as a [`data.table`][data.table::data.table-class]
#' @param in_granges `boolean`; Do you want the outcome in [`GRanges`][GenomicRanges::GRanges()] format?
#' @param order_by_sd `boolean`; Order output matrix by standard deviation
#' @param by `string`; split the matrix by `row` or `col` if `n_chunks != 1`.
#' @return If `add_loci == TRUE`, `data.frame`. If `in_granges = TRUE`, [`GRanges`][GenomicRanges::GRanges()]. Otherwise, `HDF5Matrix` or `matrix`. 
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
#' 
#' # Split the matrix into parts
#' get_matrix(scMethrix_data, n_chunks = 4, by="row")
#' @export
get_matrix <- function(scm = NULL, assay = "score", add_loci = FALSE, in_granges=FALSE, order_by_sd=FALSE, n_chunks = 1, by=c("row","column")) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  .validateType(add_loci,"boolean")
  .validateType(in_granges,"boolean")
  .validateType(order_by_sd,"boolean")
  .validateType(n_chunks,"integer")
  by <- .validateArg(by,get_matrix)
  
  if (by=="row") {
    n_chunks <- min(max(n_chunks,1),nrow(scm))
  } else {
    n_chunks <- min(max(n_chunks,1),ncol(scm))
  }
  
  if (add_loci == FALSE & in_granges == TRUE)
    warning("Without genomic locations (add_loci= FALSE), it is not possible to convert the results to GRanges, ", 
            "the output will be a data.frame object. ")
  
  if ((in_granges || add_loci) && n_chunks != 1) stop("Unable to split matrix if either in_granges = T or add_loci=T")
  
  #---- Function code ------------------------------------------------------
  mtx <- SummarizedExperiment::assay(x = scm, i = which(assay == SummarizedExperiment::assayNames(scm)))
  
  if (order_by_sd) {
    sds = DelayedMatrixStats::rowSds(mtx, na.rm = TRUE)
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
  
  if (n_chunks != 1) {
    if (by == "row") {
      row_idx <- split_vector(1:nrow(mtx),chunks=n_chunks)
      mtx <- lapply(row_idx,function(row_idx) mtx[row_idx,,drop=FALSE])
    } else {
      col_idx <- split_vector(1:ncol(mtx),chunks=n_chunks)
      mtx <- lapply(col_idx,function(col_idx) mtx[,col_idx,drop=FALSE])
    }
  }
  
  return (mtx)
}

#---- save_scMethrix ---------------------------------------------------------------------------------------------------
#' Saves an [`scMethrix-class`] object to disk
#' @description Saves an [`scMethrix-class`] object either to an `.rds` file (in-memory experiments), or into `.h5` and `.rds` files (HDF5 experiments), for assay data and experiment data, respectively. The experiments can be loaded via [load_scMethrix()].
#' @details HDF5 and in-memory [`scMethrix-class`] objects are stored differently:
#' * HDF5 objects have two files: `assays.h5` and `se.rds`. These files are hardcoded in the [HDF5Array::saveHDF5SummarizedExperiment] function and cannot be changed, as the load functions will only look for these files. To use these, you must specify the directory they will be stored in. That function is somewhat dangerously coded as well, as it will delete everything else in the directory when you save something. Here, a menu prompt has been added to warn the user.
#' * In-memory objects are simply stored in a `.rds` container. There is no requirement for file name or such, but will still prompt if the file already exists
#' Using the flag `replace = TRUE` will override the menus, but care must be taken to not delete important files (this happened numerous times to the authors!).
#' 
#' If `quick = TRUE` for HDF5 experiments, any operations done on assay matrices will not be realized. In other words, the assay information on the hard disk will not be changed. Non-matrix information will be updated (e.g., metadata) as well as any pending matrix operations. To use this, the experiment must have previously been saved using `quick = FALSE`.
#' 
#' @inheritParams generic_scMethrix_function
#' @param replace `boolean`; Should it overwrite the pre-existing data? Default = `FALSE`.
#' @param quick `boolean`; Flag to skip realizing of matrix operations
#' @param dest `string`; the destination folder for the output files
#' @param ... Additional parameters to pass to [HDF5Array::saveHDF5SummarizedExperiment()].
#' @importFrom SummarizedExperiment assays
#' @importFrom methods extends
#' @seealso [load_scMethrix()]
#' @examples
#' data('scMethrix_data')
#' dir <- paste0(tempdir(),"/h5")
#' scm <- convert_scMethrix(scMethrix_data, h5_dir=dir)
#' save_scMethrix(scm, dest = dir, replace = TRUE)
#' @return invisible [`scMethrix`] object, with the assay data stored in the `dest` (if HDF5 format)
#' @export
save_scMethrix <- function(scm = NULL, dest = NULL, replace = FALSE, quick = FALSE, verbose = TRUE, ...) {
  
  #---- Input validation ---------------------------------------------------
  if (!extends(class(scm),"SummarizedExperiment")) {
    stop("A valid SummarizedExperiment-derived object needs to be supplied.", call. = FALSE)
  }
  
  .validateType(quick,"boolean")
  .validateType(replace,"boolean")
  .validateType(verbose,"boolean")
  
  #---- Function code ------------------------------------------------------
  
  if (verbose) message("Saving scMethrix object to ",dest, start_time())
  
  if (is_h5(scm)) {
    if (quick) {
      if (!is.null(dest)) warning("dest is not used when quicksaving experiments. Experiment will be saved in it's original directory")
      exp <- HDF5Array::quickResaveHDF5SummarizedExperiment(x = scm, verbose=verbose) 
    } else {
      if (is.null(dest)) {
        dest = tempfile("scm_")
        warning("No dest specified. Experiment will be save to temporary directory:\n   ",dest) 
      } else if (dir.exists(dest) && !replace) {
        files <- list.files (dest,full.names = TRUE)
        
        if (length(files) != 0) {
          cat("Files are present in the target directory, including: \n")
          writeLines(paste("   ",head(files)))
          choice <- menu(c("Yes", "No"), title="Are you sure you want to delete this directory and save the HDF5 experiment?")
          
          if (choice == 2 || choice == 0) {
            message("Saving aborted. The target directory has not been affected.")
            return(invisible(NULL))
          } else {
            replace = TRUE
          }
        }
      }
      
      if (verbose) message("Saving HDF5 experiment to disk...",start_time())
      
      scm <- HDF5Array::saveHDF5SummarizedExperiment(x = scm, dir = dest, replace = replace, 
                                                     chunkdim = c(length(rowRanges(scm)),1), level = 6, verbose = verbose,...)
    } 
  } else {
    
    if (file.exists(dest) && !replace) {
      choice <- menu(c("Yes", "No"), title="Destination file already exists. Overwrite?")
      if (choice == 2 || choice == 0) {
        message("Saving aborted. The target directory has not been affected.")
        return(invisible(NULL))
      } 
    }
    
    saveRDS(scm,dest)
  }
  
  if (verbose) message("Experiment saved in ",stop_time())
  
  return(invisible(scm))
  
}

#---- load_scMethrix ---------------------------------------------------------------------------------------------------
#' Loads HDF5 [`scMethrix-class`] object
#' @description Loads an [`scMethrix-class`] object from the hard disk. Experiments can be saved via [save_scMethrix()]
#' @inheritParams generic_scMethrix_function
#' @param dest `string`; The `directory` or `file` to read in from
#' @param ... Parameters to pass to [HDF5Array::loadHDF5SummarizedExperiment]
#' @return An [`scMethrix-class`] object
#' @seealso [save_scMethrix()]
#' @examples
#' data('scMethrix_data')
#' dir <- paste0(tempdir(),"/h5")
#' scm <- convert_scMethrix(scMethrix_data, h5_dir=dir)
#' save_scMethrix(scm, dest = dir, replace = TRUE)
#' n <- load_scMethrix(dest = dir)
#' @export
load_scMethrix <- function(dest = NULL, verbose = TRUE, ...) {
  
  .validateType(verbose,"boolean")
  
  if (verbose) message("Loading scMethrix object", start_time())
  
  if (.validateType(dest,"file",throws=F)) {
    scm <- readRDS(dest)
  } else if (.validateType(dest,"dir",throws=F)) {
    scm <- HDF5Array::loadHDF5SummarizedExperiment(dir = dest, ...)
    scm <- as(scm, "scMethrix")
  } else {
    stop("Invalid type. Must be either a directory or file for 'dest'")
  }
  
  if (verbose) message("Loaded in ",stop_time())
  
  return(scm)
}

#--- convert_scMethrix ----------------------------------------------------------------------------------------------------
#' Converts an in-memory [`scMethrix-class`] to an HDF5 [`scMethrix-class`]
#' @details Takes a [`scMethrix-class`] object and returns with the same object with delayed array assay slots
#' with HDF5 backend. Might take long time!
#' @inheritParams generic_scMethrix_function
#' @param type `string`; what type of scMethrix to convert to. If `NULL`, this will convert to the opposite type, otherwise, will convert (if necessary) to the type specified
#' @return An [`scMethrix-class`] object in HDF5 format
#' @importFrom SummarizedExperiment assays
#' @examples
#' data('scMethrix_data')
#' convert_scMethrix(scMethrix_data, h5_dir=paste0(tempdir(),"/h5"))
#' @export
convert_scMethrix <- function(scm = NULL, type = c(NA,"HDF5","memory"), h5_dir = NULL, verbose = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  type <- .validateArg(type,convert_scMethrix)
  .validateType(h5_dir,c("string","null"))
  .validateType(verbose,"boolean")
  
  if (is.na(type)) {
    type <- ifelse(is_h5(scm),"memory","HDF5")
  }
  
  if (is_h5(scm) && type == "HDF5") return(scm)
  
  if (!is_h5(scm) && type == "memory") return(scm)
  
  #---- Function code ------------------------------------------------------
  
  if (type == "HDF5") {
    
    if (verbose) message("Converting in-memory scMethrix to HDF5", start_time())
    
    metadata <- metadata(scm)
    metadata <- metadata[names(metadata) == "is_h5"]
    
    scm <- scMethrix(assays = assays(scm), h5_dir = h5_dir, rowRanges = rowRanges(scm), is_h5 = TRUE, 
                     colData = colData(scm), replace = TRUE, verbose = verbose)
  } else {
    
    if (verbose) message("Converting HDF5 scMethrix to in-memory", start_time())
    
    for (name in SummarizedExperiment::assayNames(scm)) {
      SummarizedExperiment::assay(scm,name) <- as.matrix(SummarizedExperiment::assay(scm,name))   
      if (verbose) message("   Converted '",name,"' assay in ", split_time())
    }
    
    scm@metadata$is_h5 <- FALSE
  }
  
  if (verbose) message("Converted in ", stop_time())
  
  validObject(scm)
  return(scm)
}

#---- subset_scMethrix -------------------------------------------------------------------------------------------------
#' Subsets an [`scMethrix-class`] object based on `regions`, `contigs` and/or `samples`.
#' @description Takes [`scMethrix-class`] object and filters CpGs based on region, contig and/or sample. Can 
#' either subset (`include`) to or filter (`exclude`) the specified parameters.
#' @inheritParams generic_scMethrix_function
#' @param regions [`GRanges`][GenomicRanges::GRanges()]; genomic regions to subset by. Can also be a [`data.table`][data.table::data.table-class] that follows [this][cast_datatable()] format.
#' @param contigs `string`; array of chromosome names to subset by
#' @param samples `string`; array of sample names to subset by
#' @param by `string`; Subset to `include` or `exclude` the given criteria from the subset. Default = `include`.
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
#' @return An object of class [`scMethrix`]
#' @export
subset_scMethrix <- function(scm = NULL, regions = NULL, contigs = NULL, samples = NULL, by=c("include","exclude"), overlap_type=c("within", "start", "end", "any", "equal"),verbose=TRUE) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)  
  .validateType(regions,c("Granges","null"))
  .validateType(contigs,c("string","null"))
  .validateType(samples,c("string","null"))
  by <- .validateArg(by,subset_scMethrix)
  overlap_type <- .validateArg(overlap_type,subset_scMethrix)
  .validateType(verbose,"boolean")
  
  if (is.null(regions) & is.null(contigs) & is.null(samples)) 
    stop("At least 1 argument mandatory for subsetting. No subset generated")
  
  
  #---- Function code ------------------------------------------------------
  if (verbose) message("Subsetting CpG sites...",start_time())
  
  if (by == "exclude") {
    
    if (!is.null(regions)) {
      regions <- cast_granges(regions)
      if (verbose) message("   Subsetting by regions")
      reg <- IRanges::subsetByOverlaps(SummarizedExperiment::rowRanges(scm), regions, invert = TRUE, type=overlap_type, maxgap=-1L, minoverlap=0L)
      scm <- subset_scMethrix(scm,regions=IRanges::reduce(reg),by="include")
    }
    
    if (!is.null(contigs)) {
      if (verbose) message("   Subsetting by contigs")
      c <- as.character(GenomeInfoDb::seqnames(scm)@values)
      scm <- subset_scMethrix(scm,contigs = c[!c %in% contigs],by="include")
    }
    
    if (!is.null(samples)) {
      if (verbose) message("   Subsetting by samples")
      s <- row.names(SummarizedExperiment::colData(scm))
      scm <- subset_scMethrix(scm,samples = s[!s %in% samples],by="include")
    }
    
  } else {
    
    if (!is.null(regions)) {
      regions <- cast_granges(regions)
      if (verbose) message("   Subsetting by regions")
      scm <- scm[GenomicRanges::findOverlaps(rowRanges(scm), regions)@from]
    }
    
    if (!is.null(contigs)) {
      if (verbose) message("   Subsetting by contigs")
      scm <- subset(scm, subset = as.vector(GenomeInfoDb::seqnames(scm)) %in% contigs)
    }
    
    if (!is.null(samples)) {
      if (verbose) message("   Subsetting by samples")
      scm <- subset(scm, select = row.names(SummarizedExperiment::colData(scm)) %in% samples)
    }
    
  }
  
  if (nrow(scm) == 0) stop("Subsetting resulted in zero entries", call. = FALSE)
  
  if (verbose) message("Subset in ",stop_time())
  
  validObject(scm)
  return(scm)
  
}

#---- getStats ------------------------------------------------------------------------------------------------------------
#' Estimate descriptive statistics for each sample
#' @details Calculate descriptive statistics (`Mean`, `Median`, `SD`, `Count`) for samples or chromosomes.
#' @inheritParams generic_scMethrix_function
#' @param perSample `boolean`; Estimate stats per sample Default = `TRUE`
#' @param perChr `boolean`; Estimate stats per chromosome. Default = `TRUE`
#' @param stats `string`; the stats to include. Default is `c("Mean","Median","SD","Count")`).
#' @examples
#' data('scMethrix_data')
#' 
#' #Get stats for each sample and chromosome
#' getStats(scMethrix_data)
#' 
#' #Get stats for each sample
#' getStats(scMethrix_data, perChr = FALSE)
#' 
##' #Get stats for each chromosome
#' getStats(scMethrix_data, perSample = FALSE)
#' @return data.table of summary stats
#' @export
getStats <- function(scm = NULL, assay="score", regions = NULL, stats = c("Mean","Median","SD","Count","Sparsity"), 
                     perSample = TRUE, perChr = TRUE, phenotype = NULL, verbose = TRUE) {

  #---- Input validation ---------------------------------------------------
  .validateExp(scm)  
  assay <- .validateAssay(scm,assay)
  .validateType(perSample,"boolean")
  .validateType(perChr,"boolean")
  .validateType(verbose,"boolean")
  .validateColData(scm, phenotype = phenotype)
  stats <- .validateArg(stats, getStats, multiple.match=T)

  Chromosome <- Sample <- Count <- Mean <- SD <- Sparsity <- ..cols <- x <- . <- NULL
  
  calc_mean <- "Mean" %in% stats
  calc_median <- "Median" %in% stats
  calc_SD <- "SD" %in% stats
  calc_count <- "Count" %in% stats
  calc_sparsity <-  "Sparsity" %in% stats
  cols <- as.character(stats)
  
  #---- Function code ------------------------------------------------------
  if (verbose) message("Getting descriptive statistics...",start_time())

  mtx <- get_matrix(scm, assay = assay)
  chrs <- as.data.table(.getGRchrStats(rowRanges(scm)))
  
  stats <- lapply(1:nrow(chrs), function(x) {
    rows <- chrs$Start.idx[x]:chrs$End.idx[x]

    data.table::data.table(
      Sample = colnames(scm),
      Chromosome = chrs$Chromosome[x],
      Mean =   if (calc_mean)   DelayedMatrixStats::colMeans2(mtx, rows = rows, na.rm = TRUE),
      SD =     if (calc_SD)     DelayedMatrixStats::colSds(mtx, rows = rows, na.rm = TRUE),
      Median = if (calc_median) DelayedMatrixStats::colMedians(mtx, rows = rows, na.rm = TRUE),
      Count =  length(rows) - DelayedMatrixStats::colCounts(mtx, rows = rows, value=NA)
    )
  })
  
  stats <- data.table::rbindlist(l = stats, use.names = TRUE)

  stats <<- stats
  
  if (!perChr) {
    
    count = unlist(stats[, .(Count = lapply(.(Count),sum)), by=.(Sample)]$Count)
    
    if (calc_mean)   {
      stats[, Mean := Mean*Count]
      mean = unlist(stats[, .(Mean = lapply(.(Mean),sum)), by=.(Sample)]$Mean)/count
    }
    
    if (calc_SD)    {
      stats[, SD := (SD*SD)*(Count)]
      sd = (unlist(stats[, .(SD = lapply(.(SD),sum)), by=.(Sample)]$SD)/count)^.5
    }
    
    if (calc_median) {
      median = stats[, .SD[which.max(Count)], by=Sample]$Median
    }

    stats <- data.table::data.table(
      Sample = levels(factor(stats$Sample)),
      Chromosome = "All",
      Mean =   if (calc_mean)   mean,
      SD =     if (calc_SD)     sd,
      Median = if (calc_median) median,
      Count =  count
    )
  }
  
  if (!perSample) {
    
    count = unlist(stats[, .(Count = lapply(.(Count),sum)), by=.(Chromosome)]$Count)
    
    if (calc_mean)   {
      stats[, Mean := Mean*Count]
      mean = unlist(stats[, .(Mean = lapply(.(Mean),sum)), by=.(Chromosome)]$Mean)/count
    }
    
    if (calc_SD)    {
      stats[, SD := (SD*SD)*(Count)]
      sd = (unlist(stats[, .(SD = lapply(.(SD),sum)), by=.(Chromosome)]$SD)/count)^.5
    }
    
    if (calc_median) {
      median = stats[, .SD[which.max(Count)], by=Chromosome]$Median
    }
    
    stats <- data.table::data.table(
      Sample = "All",
      Chromosome = levels(factor(stats$Chromosome)),
      Mean =   if (calc_mean)   mean,
      SD =     if (calc_SD)     sd,
      Median = if (calc_median) median,
      Count =  count
    )
  } 

  if (calc_sparsity) {
    if (perSample && perChr) {
      chrs <- chrs[,c("Chromosome","Sites")]
      setnames(chrs, "Sites", "Sparsity")
      stats <- merge(stats,chrs,by="Chromosome")
      stats[, Sparsity := 1-(Count/Sparsity)]
    } else if (perChr) {
      props <- lapply(1:nrow(chrs), function(x) {
        rows <- chrs$Start.idx[x]:chrs$End.idx[x]
        data.table::data.table(
          Chromosome = chrs$Chromosome[x],
          Sparsity = 1-((length(rows) - sum(DelayedMatrixStats::rowAlls(mtx, rows = rows, value=NA)))/length(rows))
        )
      })
      props <- data.table::rbindlist(l = props, use.names = TRUE)
      stats <- merge(stats,props,by="Chromosome")
    } else if (perSample) {
      stats[, Sparsity := 1-(Count/nrow(scm))]
    } else {
      stats[, Sparsity := 1-(Count/prod(dim(scm)))]
    }
  }
  
  cols <- c("Sample","Chromosome",cols)
  stats <- stats[, cols, with = FALSE]

  gc()
  if (verbose) message("Finished in ", stop_time())
  
  return(stats)
}

#---- getColDataStats -------------------------------------------------------------------------------------
#' Adds descriptive statistics to colData columns in an [`scMethrix-class`] object.
#' @details Adds the mean, SD, and sample count for each sample in an [`scMethrix-class`] object. This can be accessed using `colData()`. Columns with the names of `mean`, `sd`, and `cpg` will be automatically overwritten, but `suffix` can be used to keep multiple stats columns.
#' 
#' This data will not be updated automatically for any subset, merge, bin, etc functions.
#' 
#' @inheritParams generic_scMethrix_function
#' @param suffix `string`; a suffix to add to the string
#' @param stats `list(string)`; the stats to include. Default = `c("Mean","SD","CpGs","Sparsity")`.
#' @return An [`scMethrix-class`] object
#' @examples
#' data('scMethrix_data')
#' getColDataStats(scMethrix_data)
#' @export
getColDataStats <- function(scm, assay = "score", suffix="", stats = c("Mean","SD","CpGs","Sparsity")) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  .validateAssay(scm,assay)
  .validateType(suffix,"string")
  stats <- .validateArg(stats, getColDataStats,multiple.match=T)
  
  calc_mean <- "Mean" %in% stats
  calc_sd <- "SD" %in% stats
  calc_cpgs <- "CpGs" %in% stats
  calc_sparsity <- "Sparsity" %in% stats
  
  #---- Function code ------------------------------------------------------
  
  cpgs <- nrow(scm)-DelayedMatrixStats::colCounts(get_matrix(scm = scm,assay = assay), value = NA)
  
  stats <-
    data.table::data.table(
      Mean = if (calc_mean) DelayedMatrixStats::colMeans2(get_matrix(scm = scm,assay = assay), na.rm = TRUE),
      #median = DelayedMatrixStats::colMedians(get_matrix(scm = scm,assay = assay), na.rm = TRUE),
      SD = if (calc_sd) DelayedMatrixStats::colSds(get_matrix(scm = scm,assay = assay), na.rm = TRUE),
      CpGs = if (calc_cpgs) cpgs,
      Sparsity = if (calc_sparsity) cpgs/nrow(scm)
    )
  
  #stats <- round(stats,2)
  
  colnames(stats) <- paste0(colnames(stats),suffix)
  colData <- colData(scm)[,!(colnames(colData(scm)) %in% colnames(stats)), drop=FALSE]
  colData(scm) <- cbind(colData,stats)
  
  validObject(scm)
  return(scm)
}

#---- getRowDataStats --------------------------------------------------------------------------------------------------
#' Adds descriptive statistics to metadata columns in an [`scMethrix-class`] object.
#' @details Adds the mean, median, SD, and sample count and coverage (if present) for the `rowData()` in an [`scMethrix-class`] object. This can be accessed using `mcols()`.
#' 
#' This data will not be updated automatically for any subset, merge, bin, etc functions.
#' 
#' @inheritParams generic_scMethrix_function
#' @inheritParams getColDataStats
#' @param stats list of strings; the stats to include. Default is `c("Mean","SD","Cells","Sparsity")`
#' @return An [`scMethrix-class`] object
#' @examples
#' data('scMethrix_data')
#' getRowDataStats(scMethrix_data)
#' @export
getRowDataStats <- function(scm, assay = "score", suffix="", stats = c("Mean","SD","Cells","Sparsity")) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)  
  .validateAssay(scm,assay)
  .validateType(suffix,"string")
  
  calc_mean <- "Mean" %in% stats
  calc_sd <- "SD" %in% stats
  calc_cells <- "Cells" %in% stats
  calc_sparsity <- "Sparsity" %in% stats
  #TODO: add median
  #---- Function code ------------------------------------------------------
  
  cells <- ncol(scm)-DelayedMatrixStats::rowCounts(get_matrix(scm = scm,assay = assay), value = NA)
  
  stats <-
    data.table::data.table(
      Mean = if (calc_mean) DelayedMatrixStats::rowMeans2(get_matrix(scm = scm,assay = assay), na.rm = TRUE),
      #median_meth = DelayedMatrixStats::rowMedians(get_matrix(scm = scm,assay = assay), na.rm = TRUE),
      SD = if (calc_sd) DelayedMatrixStats::rowSds(get_matrix(scm = scm,assay = assay), na.rm = TRUE),
      Cells = if (calc_cells) cells,
      Sparsity = if (calc_sparsity) (cells/ncol(scm))
    )
  
  #stats <- round(stats,2)
  
  # Set SD to zero for rows with only one CpG (as NA rows and rows with one value will give zero SD), and set rows with zero cells to NA
  if (calc_sd) {
    stats[is.na(get("SD")), ("SD") := 0]
    stats[cells == 0, ("SD") := NA]
  }
  
  colnames(stats) <- paste0(colnames(stats),suffix)
  rowData <- rowData(scm)[,!(colnames(rowData(scm)) %in% colnames(stats)), drop=FALSE]
  rowData(scm) <- cbind(rowData,stats)
  
  validObject(scm)
  return(scm)
}


#---- getRegionStats ---------------------------------------------------------------------------------------------------
#' Extracts and summarizes methylation or coverage info by regions of interest
#' @details Summarizes regions and/or groups for descriptive statistics.
#' @inheritParams generic_scMethrix_function
#' @param by `function`; mathematical function by which regions should be summarized. Can be one of the following: `mean`, `median`, `maximum`, `minimum`, `sum`, or `sd`. Default = `mean`
#' @param group `string`; a column name from sample annotation that defines groups. In this case, the number of `min_samples` will be tested group-wise.
#' @importFrom methods setClass
#' @return table of summary statistic for the given region
#' @examples
#' data('scMethrix_data')
#' 
#' # Determine global methylation of groups
#' # getRegionStats(scMethrix_data, assay = 'score', group = "Group",  by = 'mean') 
#' #TODO: update the scMethrix_data for group info
#' 
#' # Determine methylation status of chromosomes
#' #getRegionStats(scMethrix_data, assay = "score", regions = range(rowRanges(scMethrix_data)))
#' 
#' # Determine methylation status of chromosomes by groups
#' #getRegionStats(scMethrix_data, assay = "score", group = "Group", 
#' #regions = range(rowRanges(scMethrix_data)))
#' @export
getRegionStats = function (scm = NULL, assay="score", regions = NULL, group = NULL, n_chunks=1, 
                           n_threads = 1, by = c('mean', 'median', 'maximum', 'minimum', 'sum', 'sd'), 
                           overlap_type = c("within", "start", "end", "any", "equal"), verbose = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm) 
  assay <- .validateAssay(scm,assay)
  .validateType(regions,c("Granges","null"))
  .validateType(group,c("string","null"))
  .validateType(n_chunks,"integer")
  n_threads <- .validateThreads(n_threads)
  by <- .validateArg(by,getRegionStats)
  overlap_type <- .validateArg(overlap_type,getRegionStats)
  .validateType(verbose,"boolean")
  
  if (!is.null(group) && !(group %in% colnames(scm@colData))){
    stop(paste("The column name ", group, " can't be found in colData. Please provid a valid group column."), call. = FALSE)
  }
  
  if (n_chunks > nrow(scm)) {
    n_chucks <- nrow(scm)
    warning("n_chunks exceeds number of files. Defaulting to n_chunks = ",n_chunks)
  }
  
  yid  <- NULL
  
  #---- Function code ------------------------------------------------------
  if(verbose) message("Generating region summary...",start_time())
  
  
  if (!is.null(regions)) {
    regions = cast_granges(regions)
  } else { # If no region is specifed, use entire chromosomes
    regions = range(SummarizedExperiment::rowRanges(scm))
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
    if (assay == "score") {
      dat = get_matrix(scm = scm[overlap_indices$xid,], assay = "score", add_loci = TRUE)
    } else if (assay == "counts") {
      dat = get_matrix(scm = scm[overlap_indices$xid,], assay = "counts", add_loci = TRUE)
    }
  } else {
    
    #stop("Chunking not enabled")
    
    if(nrow(overlap_indices) < n_chunks){
      n_chunks <- nrow(overlap_indices)
      warning("Fewer overlaps indicies than n_chunks. Defaulting to n_chunks = ",n_chunks)
    }
    
    # if (n_chunks < n_threads) {
    #   n_threads <- n_chunks
    #   warning("n_threads < n_chunks. Defaulting to n_threads = ",n_threads)
    # }
    
    cl <- parallel::makeCluster(n_threads)
    doParallel::registerDoParallel(cl)
    
    parallel::clusterEvalQ(cl, c(library(data.table), library(scMethrix), sink(paste0("D:/Git/scMethrix/", Sys.getpid(), ".txt"))))
    parallel::clusterEvalQ(cl, expr={
      scMethrix <- setClass(Class = "scMethrix", contains = "SingleCellExperiment")
    })
    #parallel::clusterExport(cl,list('scm','scMethrix','type','is_h5','get_matrix','start_time','split_time','stop_time'))
    
    chunk_overlaps <- split(overlap_indices$xid, ceiling(seq_along(overlap_indices$xid) /
                                                           ceiling(length(overlap_indices$xid)/n_chunks)))
    
    data <- parallel::parLapply(cl, chunk_overlaps, fun = function(i) {
      get_matrix(scm[i,], assay = assay, add_loci = TRUE) # TODO: object of type 'S4' is not subsettable
    })
    
    parallel::stopCluster(cl)
    
    data <- rbindlist(data)
    
  }
  
  if(nrow(overlap_indices) != nrow(dat)){
    stop("Something went wrong")
  }
  
  dat = cbind(overlap_indices, dat)
  
  #message("-Summarizing overlaps..\n")
  if(by == "mean") {
    message("Summarizing by average")
    output = dat[, lapply(.SD, mean, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(scm)))]
  } else if (by == "median") {
    message("Summarizing by median")
    output = dat[, lapply(.SD, median, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(scm)))]
  } else if (by == "maximum") {
    message("Summarizing by maximum")
    output = suppressWarnings(
      dat[, lapply(.SD, max, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(scm)))])
    for (j in 1:ncol(output)) set(output, which(is.infinite(output[[j]])), j, NA)
  } else if (by == "mininum") {
    message("Summarizing by minimum")
    output = suppressWarnings(
      dat[, lapply(.SD, min, na.rm = TRUE), by = yid, .SDcols = c(rownames(colData(scm)))])
    for (j in 1:ncol(output)) set(output, which(is.infinite(output[[j]])), j, NA)
  } else if (by == "sum") {
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

# #--- expand_scMethrix -----------------------------------------------------------------------------------------------------
# #' Expands an [`scMethrix`] object to match input \code{regions}.
# #' @details Takes [`scMethrix`] object and adds CpGs to the object
# #' @inheritParams generic_scMethrix_function
# #' @param regions genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object
# #' @param overlap_type string; defines the type of the overlap of the CpG sites with the target region. Default value is `within`. For detailed description, see the \code{findOverlaps} function of the \code{\link{IRanges}} package.
# #' @examples
# #' data('scMethrix_data')
# #' @return An object of class [`scMethrix`]
# #' @export
# expand_scMethrix <- function(scm = NULL, regions = NULL, overlap_type=c("within", "start", "end", "any", "equal"),verbose=TRUE) {
# 
#   #---- Input validation ---------------------------------------------------
#   .validateExp(scm)
#   .validateType(regions,c("Granges","null"))
#   overlap_type <- .validateArg(overlap_type,expand_scMethrix)
#   .validateType(verbose,"boolean")
# 
#   #---- Function code ------------------------------------------------------
#   if (verbose) message("Expanding CpG sites...",start_time())
# 
#   if (!is.null(regions)) {
#     regions <- cast_granges(regions)
#     if (verbose) message("   Subsetting by regions")
#     scm <- scm[GenomicRanges::findOverlaps(rowRanges(scm), regions)@from]
#   }
# 
#   if (verbose) message("Subset in ",stop_time())
# 
#   return(scm)
# 
# }

#--- remove_uncovered ---------------------------------------------------------------------------------------
#' Remove loci that are uncovered across all samples
#' @details Takes [`scMethrix-class`] object and removes loci that are uncovered across all samples
#' @inheritParams generic_scMethrix_function
#' @return An object of class [`scMethrix-class`]
#' @examples
#' data('scMethrix_data')
#' # Remove uncovered CpGs after subsetting to a single sample
#' remove_uncovered(subset_scMethrix(scMethrix_data, samples = "C1", by="include"))
#' @export
remove_uncovered <- function(scm = NULL, n_threads = 1, verbose = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  .validateType(verbose,"boolean")
  n_threads <- .validateThreads(n_threads)
  
  #---- Function code ------------------------------------------------------
  if (verbose) message("Removing uncovered CpGs...", start_time())
  
  check_uncovered <- function (x) rowSums(!is.na(x))==0
  
  #check_uncovered <- function(mtx) DelayedMatrixStats::rowAlls(mtx,value=NA)
  
  if (n_threads != 1) {
    
    cl <- parallel::makeCluster(n_threads)
    doParallel::registerDoParallel(cl)
    
    parallel::clusterEvalQ(cl, c(library(DelayedMatrixStats)))
    #parallel::clusterExport(cl,list('read_bed_by_index','start_time','split_time','stop_time','get_sample_name'))
    
    mtx <- get_matrix(scm, n_chunks = n_threads, by="row")
    
    row_idx <- unlist(parallel::parLapply(cl,mtx,fun=check_uncovered))
    
    parallel::stopCluster(cl)
    
  } else {
    row_idx <- check_uncovered(get_matrix(scm))
  }
  
  if (sum(row_idx) == nrow(scm)) stop("All CpGs were masked. No methylation data must be present.")
  
  if (verbose) message(paste0("Removed ", format(sum(row_idx), big.mark = ","),
                              " [", round(sum(row_idx)/nrow(scm) * 100, digits = 2), "%] uncovered loci of ",
                              format(nrow(scm), big.mark = ","), " sites (",stop_time(),")"))
  
  if (!sum(row_idx) == 0) scm <- scm[!row_idx, ]
  
  validObject(scm)
  return(scm)
}

#--- mask_scMethrix -------------------------------------------------------------------------------------
#' Masks rows or columns based on some descriptive statistic
#' @details Takes [`scMethrix-class`] object and masks CpG sites based on row statistics. The sites will remain in the object and all assays will be masked. These sites can later be removed with [remove_uncovered()]. 
#'  
#'  
#'  ## Types of functions
#'  This is a very flexible function, and can do operations like:
#'  * Remove samples with few CpG sites
#'  * Remove CpGs with low coverage or are not present in many cells
#'  * Remove low variance CpG sites
#'  
#'  ## Notes
#'  * For `stat = "variance"`, a CpG that is either hypo- or hyper-methylated in all samples will have a variability of `0`, whereas a CpG that is exactly half of each will have a value of `1`. 
#'  * For `stat = "variance"` and `stat = "sd"`, CpGs present in only one sample will automatically have a `variance`/`SD` of `0`
#'  @family masking
#' @family quality control
#' @inheritParams generic_scMethrix_function
#' @param threshold numeric; The threshold value to compute from
#' @param by string; Calculate over rows or columns
#' @param stat string; The calculation to perform on each row
#' @param op string; The operator to compare the calculation to the threshold
#' @return An object of class [`scMethrix-class`]
#' @importFrom SummarizedExperiment assays assays<-
#' @seealso [mask_by_idx()], the actual masking function
#' @examples
#' data('scMethrix_data')
#' 
#' ## Mask by low CpG count
#' # This will mask any samples that have < 170 total CpGs present
#' mask_scMethrix(scMethrix_data, assay="score", threshold = 170, by = "col", stat="count", op="<")
#' 
#' #' ## Mask by low sample count
#' # This will remove CpGs that are not present in > 2 samples
#' mask_scMethrix(scMethrix_data, assay="score", threshold=2, by = "row", stat="count", op="<=")
#' 
#' ## Mask by low coverage
#' # This will mask each CpG in the counts assay where the total coverage is < 5
#' mask_scMethrix(scMethrix_data, assay="counts", threshold=5, by = "row", stat="sum", op="<")
#' 
#' ## Mask by high average coverage
#' # This will mask each CpG in the counts assay where the average coverage is > 1
#' mask_scMethrix(scMethrix_data, assay="counts", threshold=1, by = "row", stat="mean", op=">")
#' 
#' ## Mask by variance
#' # This will mask low variance (homogenous) CpG sites in the score assay where var < 0.05
#' mask_scMethrix(scMethrix_data, assay="score", threshold=0.05, by = "row", stat="var", op="<")
#' @export
mask_scMethrix <- function(scm = NULL, assay="score", threshold = 0, by=c("row","column"), stat = c("count","sum","mean","median","sd","variance","proportion"), op = c("==","!=",">","<",">=","<="), na.rm = TRUE, n_threads=1 , verbose = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  stat = .validateArg(stat,mask_scMethrix)
  op = .validateArg(op,mask_scMethrix,partial.match = F)
  by = .validateArg(by,mask_scMethrix)
  .validateType(na.rm,"boolean")
  n_threads <- .validateThreads(n_threads)
  .validateType(verbose,"boolean")
  .validateType(threshold,"numeric")
  
  if (!is_h5(scm) && n_threads != 1) 
    stop("Parallel processing not supported for a non-HDF5 scMethrix object due to probable high memory usage.
         \nNumber of cores (n_threads) needs to be 1.", call. = FALSE)
  
  row_idx <- col_idx <- mtx <- NULL
  
  #---- Function code ------------------------------------------------------
  if (verbose) message("Masking ",by,"s in the '",assay,"' assay by ",stat," ",op," ",threshold,start_time())
  
  if (by == "row") {
    if (stat == "sum")        {calc <- function (mtx) DelayedMatrixStats::rowSums2  (mtx, na.rm = na.rm)}
    else if (stat == "mean")  {calc <- function (mtx) DelayedMatrixStats::rowMeans2 (mtx, na.rm = na.rm)}
    else if (stat == "mode")  {calc <- function (mtx) DelayedMatrixStats::rowMedians(mtx, na.rm = na.rm)}
    else if (stat == "count") {calc <- function (mtx) {
      ncol(mtx) - DelayedMatrixStats::rowCounts (mtx, na.rm = na.rm, value = as.integer(NA))}
    } else if (stat == "sd")    {calc <- function (mtx) {
      vals <- DelayedMatrixStats::rowSds(mtx, na.rm = na.rm)
      vals[which(is.na(vals))] = 0 # Since CpGs rep'd by a single sample return NA instead of 0
      return(vals)
    }
    } else if (stat == "variance") {calc <- function (mtx) {
      vals <- DelayedMatrixStats::rowVars(mtx, na.rm = na.rm)
      vals[which(is.na(vals))] = 0 # Since CpGs rep'd by a single sample return NA instead of 0
      return(vals*2) # Since the variance can be at most 0.5 due to beta ratio
    }
    } else if (stat =="proportion") {calc <- function(mtx) {
      (ncol(mtx) - DelayedMatrixStats::rowCounts(mtx, na.rm = na.rm, value = as.integer(NA))) / ncol(mtx)
    }
    } else {stop("Unknown calculation. This should not happen", call. = FALSE)}
    
    expr <- call(op,calc(get_matrix(scm = scm, assay=assay)),threshold)
    row_idx <- which(eval(expr))
  } else {
    if (stat == "sum")        {calc <- function (mtx) DelayedMatrixStats::colSums2(  mtx, na.rm = na.rm)}
    else if (stat == "mean")  {calc <- function (mtx) DelayedMatrixStats::colMeans2( mtx, na.rm = na.rm)}
    else if (stat == "mode")  {calc <- function (mtx) DelayedMatrixStats::colMedians(mtx, na.rm = na.rm)}
    else if (stat == "count") {calc <- function (mtx) {
      nrow(mtx) - DelayedMatrixStats::colCounts (mtx, na.rm = na.rm, value = as.integer(NA))}
    } else if (stat == "sd")    {calc <- function (mtx) {
      vals <- DelayedMatrixStats::rowSds(mtx, na.rm = na.rm)
      vals[which(is.na(vals))] = 0 # Since CpGs rep'd by a single sample return NA instead of 0
      return(vals)
    }
    } else if (stat == "variance") {calc <- function (mtx) {
      vals <- DelayedMatrixStats::colVars(mtx, na.rm = na.rm)
      vals[which(is.na(vals))] = 0 # Since CpGs rep'd by a single sample return NA instead of 0
      return(vals*2) # Since the variance can be at most 0.5 due to beta ratio
    }
    } else if (stat =="proportion") {calc <- function(mtx) {
      (nrow(mtx) - DelayedMatrixStats::colCounts(mtx, na.rm = na.rm, value = as.integer(NA))) / nrow(mtx)}
    } else {stop("Unknown calculation. This should not happen", call. = FALSE)}
    
    expr <- call(op,calc(get_matrix(scm = scm, assay=assay)),threshold)
    col_idx <- which(eval(expr))
  }
  
  scm <- mask_by_idx(scm, col_idx, row_idx, verbose = verbose)
  
  validObject(scm)
  return(scm)
}

#' Specified rows and columns are masked with NA values
#' @details This iterates through all assays in the inputted [`scMethrix-class`] object and replaces all rows in row_idx and cols in col_idx. 
#' @family masking
#' @family quality control
#' @inheritParams generic_scMethrix_function
#' @param col_idx numeric; A vector of column indexes to replace all values with NA
#' @param row_idx numeric; A vector of row indexes to replace all values with NA
#' @return An object of class [`scMethrix-class`]
#' @importFrom SummarizedExperiment assays assays<-
#' @seealso [mask_scMethrix()] wraps this function and generates idxs by statistics, [remove_uncovered()] to remove the masked sites
#' @export
mask_by_idx <- function (scm, col_idx = NULL, row_idx = NULL, verbose = TRUE) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  #.validateType(row_idx,"integer")
  .validateType(verbose,"boolean")
  
  if (length(col_idx) == ncol(scm)) stop("No samples left after masking.", call. = FALSE)
  if (length(row_idx) == nrow(scm)) stop("No CpG sites left after masking.", call. = FALSE)
  
  if (!is.null(row_idx)) {
    for (i in 1:length(assays(scm))) {
      assays(scm)[[i]][row_idx,] <- as.integer(NA)
    } 
    if (verbose) message("   Masked ",length(row_idx)," [",round(length(row_idx)/nrow(scm)*100, digits = 2), "%] CpG sites")
  }
  
  if (!is.null(col_idx)) {
    for (i in 1:length(assays(scm))) {
      assays(scm)[[i]][,col_idx] <- as.integer(NA)
    } 
    if (verbose) message("   Masked ",length(col_idx)," [",round(length(col_idx)/ncol(scm)*100, digits = 2), "%] samples")
  }
  
  if (verbose) message("Masking completed in ",stop_time())
  
  return(scm)
}
