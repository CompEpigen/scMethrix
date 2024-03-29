% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scmethrix_operations.R
\name{remove_assay}
\alias{remove_assay}
\title{Removes an assay from an \code{\link{scMethrix}} object}
\usage{
remove_assay(scm = NULL, assay = NULL)
}
\arguments{
\item{scm}{\code{\link{scMethrix}}; the single cell methylation experiment}

\item{assay}{string; name of an existing assay. Default = "score"}
}
\value{
An \code{\link{scMethrix}} object
}
\description{
Removes an assay from an \code{\link{scMethrix}} object
}
\details{
This will remove an assay from the scMethrix experiment object. All transformed assays may be removed, as well as the coverage assay (since it is less useful when compared to normal WGBS data), but the score assay cannot be removed. Reduced dimensionality data will be retained even if the parent assay is removed.
}
\examples{
data('scMethrix_data')
remove_assay(scMethrix_data,assay="counts")
}
