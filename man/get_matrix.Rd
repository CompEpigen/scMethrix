% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scmethrix_operations.R
\name{get_matrix}
\alias{get_matrix}
\title{Extract assays from an \code{\link{scMethrix}} object}
\usage{
get_matrix(
  scm = NULL,
  assay = "score",
  add_loci = FALSE,
  in_granges = FALSE,
  order_by_sd = FALSE,
  n_chunks = 1,
  by = c("row", "column")
)
}
\arguments{
\item{scm}{\code{\link{scMethrix}}; the single cell methylation experiment}

\item{assay}{string; name of an existing assay. Default = "score"}

\item{add_loci}{Default FALSE. If TRUE adds CpG position info to the matrix and returns as a data.table}

\item{in_granges}{Do you want the outcome in \code{\link{GRanges}}?}

\item{order_by_sd}{Order output matrix by standard deviation}

\item{n_chunks}{integer; Number of chunks to split the \code{\link{scMethrix}} object in case it is very large. Default = 1}

\item{by}{string; split the matrix by "row" or "col" if n_chunks != 1}
}
\value{
HDF5Matrix or matrix
}
\description{
Extract assays from an \code{\link{scMethrix}} object
}
\details{
Takes \code{\link{scMethrix}} object and returns the \code{methylation} matrix. This will return in the format used by the object (matrix or HDF5matrix).
}
\examples{
data('scMethrix_data')

# Get methylation data
get_matrix(scMethrix_data)

# Get methylation data with loci
get_matrix(scMethrix_data, add_loci=TRUE)

# Get methylation data with loci inside a Granges object 
get_matrix(scMethrix_data, add_loci=TRUE, in_granges=TRUE)

# Get methylation data sorted by SD
get_matrix(scMethrix_data, order_by_sd = TRUE)

# Split the matrix into parts
get_matrix(scMethrix_data, n_chunks = 4, by="row")
}
