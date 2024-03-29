% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scmethrix_operations.R
\name{mask_helper}
\alias{mask_helper}
\title{Helper function for masking. All rows in row_idx will be set to NA.}
\usage{
mask_helper(scm, row_idx, verbose = TRUE)
}
\arguments{
\item{scm}{\code{\link{scMethrix}}; the single cell methylation experiment}

\item{row_idx}{numeric; A vector of row indexes for which to replace with NA}

\item{verbose}{boolean; Flag for outputting function status messages. Default = TRUE}
}
\value{
An object of class \code{\link{scMethrix}}
}
\description{
Helper function for masking. All rows in row_idx will be set to NA.
}
\details{
Used with mask_by_sample, mask_by_coverage, and mask_by_variance. It iterates through all assays in the inputted \code{\link{scMethrix}} object and replaces all rows in row_idx with NA.
}
