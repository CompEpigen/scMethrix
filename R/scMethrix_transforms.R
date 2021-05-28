#------------------------------------------------------------------------------------------------------------
#' Transforms an assay in an \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @param scm A \code{\link{scMethrix}} object
#' @param assay String name of an existing assay
#' @param name String name of transformed assay
#' @param trans The transformation function
#' @param h5_temp temporary directory to store hdf5
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' transform_assay(scMethrix_data,assay="score",name="plus1",trans=function(x){x+1})
#' @export
transform_assay <- function(scm,assay = NULL, name = NULL, trans = NULL, h5_temp = NULL) {
  
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (typeof(trans) != "closure") {
    stop("A valid transform function must be specified.", call. = FALSE)
  }
  
  if (!(assay %in% SummarizedExperiment::assayNames(scm))) {
    stop("Assay does not exist in the object", call. = FALSE)
  }
  
  if (name %in% SummarizedExperiment::assayNames(scm)) {
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)
  }
  
  if (is_h5(scm)) {
    
    if (is.null(h5_temp)) {h5_temp <- tempdir()}
    
    grid <- DelayedArray::RegularArrayGrid(refdim = dim(scm),
                                           spacings = c(length(scm), 1L)) 
    
    trans_sink <- HDF5Array::HDF5RealizationSink(dim = dim(scm),
                                                 dimnames = list(NULL,row.names(colData(scm))), type = "integer",
                                                 filepath = tempfile(pattern="trans_sink_",tmpdir=h5_temp),
                                                 name = name, level = 6)
    
    blocs <- DelayedArray::blockApply(get_matrix(scm,type=assay), grid = grid, FUN = trans)
    
    for(i in 1:length(blocs)) {
      
      DelayedArray::write_block(block = as.matrix(blocs[[i]]), viewport = grid[[i]], sink = trans_sink)
      
    }
    
    rm(blocs)
    s <- as(trans_sink, "HDF5Array")
    
  } else {
    s <- apply(get_matrix(scm,type=assay),1:2,trans)
  }
  
  assays(scm)[[name]] <- s
  
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Bins the ranges of an \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @param scm A \code{\link{scMethrix}} object
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
bin_scMethrix <- function(scm, bin_size = 100000, trans = NULL, h5_dir = NULL) {

  if (is.null(trans)) {
    trans <- c(score = function(x) mean(x,na.rm=TRUE),counts = function(x) sum(x,na.rm=TRUE))
  }
    
  if (!is(scm, "scMethrix")) {
    stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
  }
  
  if (is_h5(scm) && is.null(h5_dir)) stop("Output directory must be specified", call. = FALSE)
  
  bins <- bin_granges(rowRanges(scm),bin_size = bin_size)
  bins <- subsetByOverlaps(bins,rowRanges(scm))
  
  sites <- sapply(1:length(bins),function (i) {
    (findOverlaps(rowRanges(scm),bins[i]))@from
  })
  
  assys <- list()
  
  for (n in 1:length(assays(scm))) {

    name <- SummarizedExperiment::assayNames(scm)[n]
    
    tryCatch(
      expr = {op <- trans[[name]]},
      error = function(e){op <- function(x) mean(x,na.rm=TRUE)}
    )
    
    vals <- lapply(sites,function (i) {
      apply(get_matrix(scm[i,],type=name),2,op)
    })
    
    setDT(vals, key=names(vals[[1]]))
    vals <- data.table::transpose(vals)
    colnames(vals) <- rownames(colData(scm))
    
    assys[[name]] <- as(vals,class(get_matrix(scm,type=name)))
  }
  
  if (is_h5(scm)) {
    
    m_obj <- create_scMethrix(assays = assys, rowRanges=bins, is_hdf5 = TRUE, 
                            h5_dir = h5_dir, genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),
                            replace = replace)
  
  } else {
    
    m_obj <- create_scMethrix(assays = assys, rowRanges=bins, is_hdf5 = FALSE, 
                              genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),)
    
  }
}

#------------------------------------------------------------------------------------------------------------
#' Imputes the NA values of a \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @param scm A \code{\link{scMethrix}} object
#' @param threshold The value for cutoff in the "score" assay to determine methylated or unmethylated status. 
#' Default = 50
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' @export
impute_scMethrix <- function (scm, threshold = 50) {
  
  #impute <- transform_assay(m,assay = "score",name = "impute",trans = function(x) {ifelse(m > threshold,1,0)})
  
  
  
}




#------------------------------------------------------------------------------------------------------------
#' Generates UMAP for scMethrix
#' @details Does UMAP stuff
#' @param scm A \code{\link{scMethrix}} object
#' @return An \code{\link{scMethrix}} object
#' @import umap
#' @import ggplot2
#' @examples
#' data('scMethrix_data')
#' @export
umap_scMethrix <- function(scm) {
  
  x <- y <- NULL
  
  start_time()
  
  scm.bin <- transform_assay(scm,assay="score",name="binary",trans=binarize)
  scm.umap <- umap(t(get_matrix(scm.bin,type="binary")),n_neighbors=min(100,10))#ncol(scm)))
 
  stop_time()
  
  df <- data.frame(x = scm.umap$layout[,1],
                   y = scm.umap$layout[,2])
  
  ggplot(df, aes(x, y)) +
    geom_point()
  

}

