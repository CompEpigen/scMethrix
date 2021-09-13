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