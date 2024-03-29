% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scmethrix_operations.R
\name{subset_scMethrix}
\alias{subset_scMethrix}
\title{Subsets an \code{\link{scMethrix}} object based on \code{regions}, \code{contigs} and/or \code{samples}.}
\usage{
subset_scMethrix(
  scm = NULL,
  regions = NULL,
  contigs = NULL,
  samples = NULL,
  by = c("include", "exclude"),
  overlap_type = c("within", "start", "end", "any", "equal"),
  verbose = TRUE
)
}
\arguments{
\item{scm}{\code{\link{scMethrix}}; the single cell methylation experiment}

\item{regions}{genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a \code{GenomicRanges} object}

\item{contigs}{string; array of chromosome names to subset by}

\item{samples}{string; array of sample names to subset by}

\item{by}{string to decide whether to "include" or "exclude" the given criteria from the subset}

\item{overlap_type}{string; defines the type of the overlap of the CpG sites with the target region. Default value is \code{within}. For detailed description, see the \code{findOverlaps} function of the \code{\link{IRanges}} package.}

\item{verbose}{boolean; Flag for outputting function status messages. Default = TRUE}
}
\value{
An object of class \code{\link{scMethrix}}
}
\description{
Subsets an \code{\link{scMethrix}} object based on \code{regions}, \code{contigs} and/or \code{samples}.
}
\details{
Takes \code{\link{scMethrix}} object and filters CpGs based on region, contig and/or sample. Can
either subset (\code{include}) to or filter (\code{exclude}) the specified parameters.
}
\examples{
data('scMethrix_data')

contigs <- c("chr1","chr3")
regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,100000000)) 
samples <- c("C1","C2")

#Subset to only samples bed1 and bed3, and chromosome 1
subset_scMethrix(scMethrix_data, samples = samples, contigs = contigs, by = "include")

#Subset to only region "chr1:1-5"
subset_scMethrix(scMethrix_data, regions = regions, by = "include")

#Subset to exclude samples bed1 and bed3, and chromosome 1
subset_scMethrix(scMethrix_data, samples = samples, contigs = contigs, by = "exclude")

#Subset to exclude region "chr1:1-5"
subset_scMethrix(scMethrix_data, regions = regions, by = "exclude")
}
