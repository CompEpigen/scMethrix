#------------------------------------------------------------------------------------------------------------
#' Transforms an assay in an \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @param m A \code{\link{scMethrix}} object
#' @param assay String name of an existing assay
#' @param name String name of transformed assay
#' @param trans The transformation function
#' @param h5_temp temporary directory to store hdf5
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' transform_assay(scMethrix_data,assay="score",name="plus1",trans=function(x){x+1})
#' @export
transform_assay <- function(m,assay = NULL, name = NULL, trans = NULL, h5_temp = NULL) {
  
  if (!is(m, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (typeof(trans) != "closure") {
    stop("A valid transform function must be specified.", call. = FALSE)
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(m))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }
  
  if (name %in% SummarizedExperiment::assayNames(m)) {
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)
  }
  
  if (is_h5(m)) {
    
    if (is.null(h5_temp)) {h5_temp <- tempdir()}
    
    grid <- DelayedArray::RegularArrayGrid(refdim = dim(m),
                                           spacings = c(dim(m)[1], 1L)) 
    
    trans_sink <- HDF5Array::HDF5RealizationSink(dim = dim(m),
                                                 dimnames = list(NULL,row.names(colData(m))), #type = "double",
                                                 filepath = tempfile(pattern="trans_sink_",tmpdir=h5_temp),
                                                 name = "trans", level = 6)
    
    blocs <- DelayedArray::blockApply(get_matrix(m,type=assay), grid = grid, trans)
    
    for(i in 1:length(blocs)) {
      
      DelayedArray::write_block(block = as.matrix(blocs[[i]]), viewport = grid[[i]], sink = trans_sink)
      
    }
    
    rm(blocs)
    s <- as(trans_sink, "HDF5Array")
    
  } else {
    s <- apply(get_matrix(m,type=assay),1:2,trans)
  }
  
  assays(m)[[name]] <- s
  
  return(m)
}

#------------------------------------------------------------------------------------------------------------
#' Bins the ranges of an \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @param m A \code{\link{scMethrix}} object
#' @param bin_size The size of each bin. First bin will begin at the start position of the first genomic
#' region on the chromosome
#' @param trans The transforms for each assay. Must be a named vector of functions (closure). 
#' Default = mean(.., na.rm=TRUE)
#' @param h5_dir directory to store H5 based object
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' trans <- c(score = function(x) mean(x,na.rm=TRUE),counts = function(x) sum(x,na.rm=TRUE))
#' bin_scMethrix(scMethrix_data,trans = trans)
#' @export
bin_scMethrix <- function(m, bin_size = 100000, trans = NULL, h5_dir = NULL) {

  if (is.null(trans)) {
    trans <- c(score = function(x) mean(x,na.rm=TRUE),counts = function(x) sum(x,na.rm=TRUE))
  }
    
  if (!is(m, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is_h5(m) && is.null(h5_dir)) stop("Output directory must be specified", call. = FALSE)
  
  bins <- bin_granges(rowRanges(m),bin_size = bin_size)
  bins <- subsetByOverlaps(bins,rowRanges(m))
  
  sites <- sapply(1:length(bins),function (i) {
    (findOverlaps(rowRanges(m),bins[i]))@from
  })
  
  assys <- list()
  
  for (n in 1:length(assays(m))) {

    name <- SummarizedExperiment::assayNames(m)[n]
    
    tryCatch(
      expr = {op <- trans[[name]]},
      error = function(e){op <- function(x) mean(x,na.rm=TRUE)}
    )
    
    vals <- lapply(sites,function (i) {
      apply(get_matrix(m[i,],type=name),2,op)
    })
    
    setDT(vals, key=names(vals[[1]]))
    vals <- data.table::transpose(vals)
    colnames(vals) <- rownames(colData(m))
    
    assys[[name]] <- as(vals,class(get_matrix(m,type=name)))
  }
  
  if (is_h5(m)) {
    
    m_obj <- create_scMethrix(assays = assys, rowRanges=bins, is_hdf5 = TRUE, 
                            h5_dir = h5_dir, genome_name = m@metadata$genome,desc = m@metadata$desc,colData = colData(m),
                            replace = replace)
  
  } else {
    
    m_obj <- create_scMethrix(assays = assys, rowRanges=bins, is_hdf5 = FALSE, 
                              genome_name = m@metadata$genome,desc = m@metadata$desc,colData = colData(m),)
    
  }
}

#------------------------------------------------------------------------------------------------------------
#' Imputes the NA values of a \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @param m A \code{\link{scMethrix}} object
#' @param threshold The value for cutoff in the "score" assay to determine methylated or unmethylated status. 
#' Default = 50
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' @export
impute_scMethrix <- function (m, threshold = 50) {
  
  #impute <- transform_assay(m,assay = "score",name = "impute",trans = function(x) {ifelse(m > threshold,1,0)})
  
  
  
}


