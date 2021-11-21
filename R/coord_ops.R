#--- liftover_CpGs ---------------------------------------------------------------------------------------------------------
#' Converts from one reference to another via rtracklayer::liftOver
#' @details 
#' This conversion is one-to-many, so for consistency, only the first element in the target assembly is used.
#' 
#' LiftOver chains can be found here: \href{https://hgdownload.soe.ucsc.edu/downloads.html}{UCSC}
#' 
#' Here's links to common chains:
#' {http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz}{hg38 to hg19}
#' {https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz}{hg19 to hg38}
#' 
#' AnnotationHub can also be used to download chains. See an example in the vignettes: 
#' 
#' @inheritParams generic_scMethrix_function
#' @param chain string; the file location for the desired chain to use in liftOver
#' @param target_genome string; the target genome. This will be update the genome field in the output scMethrix object
#' @importFrom rtracklayer import.chain liftOver
#' @export
#' @return a list of data.table containing number of CpG's and contig lengths
#' @examples
#'\dontrun{
#' hg19_cpgs = methrix::extract_CPGs(ref_genome = 'BSgenome.Hsapiens.UCSC.hg19')
#' }
liftover_CpGs <- function(scm, chain = NULL, target_genome = NULL, verbose = TRUE) {
  
  #- Input Validation --------------------------------------------------------------------------
  
  .validateExp(scm)
  .validateType(target_genome,"string")
  
  if (verbose) message("Applying liftover...")
  
  #- Function code -----------------------------------------------------------------------------
  
  n_cpg = nrow(scm)
  rrng.new <- rtracklayer::liftOver(rowRanges(scm),chain)
  
  # Remove the missing probes
  row_idx <- (lengths(rrng.new) == 0)
  scm <- scm[!row_idx,]
  rrng.new <- rrng.new[!row_idx,]
  
  # Select the first coord of one-to-many coords
  row_idx <- which(lengths(rrng.new) > 1)
  rrng.new[row_idx] <- S4Vectors::endoapply(rrng.new[row_idx],head,1)
  
  if (verbose) message("Lost ", n_cpg - nrow(scm), " CpGs during liftOver." )
  
  rowRanges(scm) <- unlist(rrng.new)
  scm@metadata$genome <- target_genome
  
  return(scm)
}

#--- extract_CpGs ---------------------------------------------------------------------------------------------------------
#' Extracts all CpGs from a genome
#' @param ref_genome BSgenome object or name of the installed BSgenome package. Example: BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome installed.genomes getBSgenome seqnames
#' @export
#' @return a list of data.table containing number of CpG's and contig lengths
#' @examples
#'\dontrun{
#' hg19_cpgs = methrix::extract_CPGs(ref_genome = 'BSgenome.Hsapiens.UCSC.hg19')
#' }
extract_CpGs = function(ref_genome = NULL) {
  
  #- Input Validation --------------------------------------------------------------------------
  
  .validateType(ref_genome,"string")
  pkgname <- seqlengths <- chr <- NULL
  
  gnoms_installed = BSgenome::installed.genomes(splitNameParts = TRUE)
  data.table::setDT(x = gnoms_installed)
  
  if (is.null(ref_genome)) {
    if (nrow(gnoms_installed) == 0) {
      stop("Could not find any installed BSgenomes.\nUse BSgenome::available.genomes() for options.")
    } else {
      message("Found following BSgenome installations. Use the required 'pkgname'.")
      print(gnoms_installed)
      stop()
      # ref_genome = gnoms_installed[,pkgname][1]
    }
  } else {
    if (nrow(gnoms_installed[pkgname %in% ref_genome]) == 0) {
      message("Could not find BSgenome ", ref_genome)
      if (nrow(gnoms_installed) == 0) {
        stop("Could not find any installed BSgenomes either.\nUse BSgenome::available.genomes() for options.")
      } else {
        message("Found following BSgenome installations. Provide the correct 'pkgname'")
        print(gnoms_installed)
        stop()
      }
    }
  }
  requireNamespace(ref_genome, quietly = TRUE)
  
  #- Function code -----------------------------------------------------------------------------
  
  message("Extracting CpGs from ",ref_genome,"...",start_time())
  
  ref_genome = BSgenome::getBSgenome(genome = ref_genome)
  
  chrs = GenomeInfoDb::standardChromosomes(ref_genome)
  
  # Code borrowed from from: https://support.bioconductor.org/p/95239/
  cgs = lapply(chrs, function(x) start(Biostrings::matchPattern("CG", ref_genome[[x]])))
  cpgs = do.call(c, lapply(seq_along(chrs), function(x) GenomicRanges::GRanges(names(ref_genome)[x],
                                                                               IRanges::IRanges(cgs[[x]], width = 2))))
  cpgs = data.table::as.data.table(as.data.frame(cpgs, stringsAsFactors = FALSE))
  cpgs[,width:=NULL][,strand:=NULL]
  colnames(cpgs) = c("chr", "start", "end")
  cpgs[, `:=`(chr, as.character(chr))][, `:=`(start, as.numeric(start))][, `:=`(end,
                                                                                as.numeric(end))]
  data.table::setkey(x = cpgs, "chr", "start")
  message(paste0("Extracted ", format(nrow(cpgs), big.mark = ","), " CpGs from ",
                 length(chrs), " contigs (",stop_time(),")"))
  
  return(cpgs)
}

#--- subset_ref_cpgs ---------------------------------------------------------------------------------------------------------
#' Subsets a given list of CpGs by another list of CpGs
#' @details Typically used to reduce the number of potential CpG sites to include only those present  in the input files so as to maximize performance and minimize resources. Can also be used for quality control to see if there is excessive number of CpG sites that are not present in the reference genome.
#' @param ref_cpgs data.table; A reference set of CpG sites (e.g. Hg19 or mm10) in bedgraph format
#' @param gen_cpgs data.table; A subset of CpG sites. Usually obtained from \code{\link{read_index}}.
#' @param verbose boolean; flag to output messages or not
#' @return Returns list of CpG sites in bedgraph format
#' @examples
#' ref_cpgs = data.frame(chr="chr1",start=(1:5*2-1), end=(1:5*2))
#' subset_ref_cpgs(ref_cpgs,ref_cpgs[1:3,])
#' @export
subset_ref_cpgs <- function(ref_cpgs, gen_cpgs, verbose = TRUE) {
  
  #- data.table initialization -----------------------------------------------------------------
  id <- NULL
  
  #- Function code -----------------------------------------------------------------------------
  keys <- rbind(ref_cpgs[,c("chr","start")], gen_cpgs[,c("chr","start")])
  data.table::setDT(keys)[, id := .GRP, by = c("chr","start")]
  
  ref <- nrow(ref_cpgs)
  gen <- nrow(gen_cpgs)
  
  keys <- list(
    ref = keys[seq_len(ref),id],
    sub = keys[ref + seq_len(gen),id]
  )
  
  sub_cpgs <- ref_cpgs[keys$ref %in% keys$sub, , drop = FALSE]
  sub <- nrow(sub_cpgs)
  
  if (verbose) message("Dropped ",ref-sub,"/",ref," CpGs (",round((ref-sub)/ref*100,2),"%) from the reference set")
  if (verbose) message(gen-sub,"/",gen," subset CpGs (",round((gen-sub)/gen*100,2),"%) were not present in the reference set")
  
  return(sub_cpgs)
}

#--- bin_granges ---------------------------------------------------------------------------------------------------------
#' Bins each region in a \code{\link{GRanges}} object into bins of specified \code{bin_size} 
#' @details Bins a single region in \code{\link{GRanges}} format into multiple regions with a specified \code{bin_size}. If \code{length(gr) %% bin_size != 0}, then the last GRange will have a length < \code{bin_size}. This is used instead of tile when you need consistently sized bins with the last bin being smaller
#' @param gr GRanges; The \code{\link{GRanges}} object
#' @param bin_size integer; x > 0; The length of region in each bin
#' @return \code{\link{GRangesList}} containing all the binned \code{\link{GRanges}}
#' @import GenomicRanges
#' @examples
#' regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,10000))
#' bin_granges(regions,bin_size=1000) 
#' @export
bin_granges <- function(gr, bin_size = 100000) {#, enforce_size = FALSE) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateType(gr,"Granges")
  .validateType(bin_size, "integer")
  .validateValue(bin_size,">0")
  
  #- Function code -----------------------------------------------------------------------------
  gr <- GenomicRanges::slidingWindows(gr,width=bin_size,step=bin_size)
  return(unlist(as(gr, "GRangesList")))
}

#--- cast_granges ---------------------------------------------------------------------------------------------------------
#' Casts genomic regions into \code{\link{GRanges}} format
#' @details Casts the input as a \code{\link{GRanges}} object. Input can be \code{\link{GRanges}} or a 
#' \code{\link{data.frame}}-compatible class that can be cast through \code{as.data.frame()}. Input BED format
#'  must be \code{chr-start-end} for \code{\link{data.frame}} objects.
#' @param regions GRanges or data.frame; The input regions to cast to \code{\link{GRanges}}
#' @return \code{\link{GRanges}} object with the input regions
#' @import GenomicRanges
#' @examples
#' regions = data.table(chr = 'chr1', start = 1, end = 100)
#' cast_granges(regions) 
#' @export
cast_granges <- function(regions) {
  if (is(regions, "GRanges")) { # do nothing
  } else if (is(regions,"data.frame")) {return (GenomicRanges::makeGRangesFromDataFrame(regions))
  } else {stop("Invalid input class for regions. Must be a GRanges or data.frame-like")}
  return(regions)
}

#--- cast_datatable ---------------------------------------------------------------------------------------------------------
#' Casts genomic regions into \code{\link{data.table}} format
#' @details Casts the input as a \code{\link{data.table}} object. Input can be \code{\link{GRanges}} or a 
#' \code{\link{data.frame}}-compatible class that can be cast through \code{as.data.frame()}. Input BED format
#'  must be \code{chr-start-end} for \code{\link{data.frame}} objects.
#' @param regions GRanges or data.frame; The input regions to cast to \code{\link{GRanges}}
#' @return \code{\link{data.table}} object with the input regions
#' @import GenomicRanges
#' @examples
#' regions = data.table(chr = 'chr1', start = 1, end = 100)
#' cast_granges(regions) 
#' @export
cast_datatable <- function(regions) {
  if (is(regions,"data.frame")) { # do nothing
  } else if (is(regions, "GRanges")) {
    regions <- as.data.table(regions)
    regions[,width:=NULL]
    setnames(regions, "seqnames", "chr")
  } else {
    stop("Invalid input class for regions. Must be a GRanges or data.frame-like")}
  return(regions)
}