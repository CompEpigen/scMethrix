% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scMethrix_transforms.R
\name{collapse_samples}
\alias{collapse_samples}
\title{Collapses multiple samples into a single sample by group}
\usage{
collapse_samples(
  scm = NULL,
  colname = NULL,
  trans = NULL,
  h5_dir = NULL,
  batch_size = 1e+05,
  n_threads = 1,
  replace = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{scm}{\code{\link{scMethrix}}; the single cell methylation experiment}

\item{colname}{string; The colname from \code{colData(scm)} indicating which samples should be collapse together}

\item{trans}{named vector of closures; The transforms for each assay in a named vector. Default NULL, meaning that
operations for "counts" assay is sum(x, na.rm=TRUE), and for all other assays is mean(x, na.rm=TRUE)}

\item{h5_dir}{string; The directory to use. Will be created if it does not exist. Default = NULL}

\item{batch_size}{The number of CpGs to calculate at once.
\code{\link{IRanges}} package.}

\item{n_threads}{integer; Maximum number of parallel instances. Default = 1}

\item{replace}{boolean; flag for whether to delete the contents of h5_dir before saving}

\item{verbose}{boolean; Flag for outputting function status messages. Default = TRUE}
}
\value{
An \code{\link{scMethrix}} object
}
\description{
Collapses multiple samples into a single sample by group
}
\details{
Multiple samples can be collapsed into a single meta-sample. Grouping for samples can be defined via colData. The collapse function can accept an arbitrary function for each assay on how to handle the collapsing (typically \code{mean} for scores, and \code{sum} for counts).

In the output object, \code{colData()} will contain a comma-delimited list of samples (\code{Samples}) that each group contains as well as the total number of CpGs in the group (\code{n_Samples}).

Reduced dimensionality data will be discarded.
}
\examples{
data('scMethrix_data')
colData(scMethrix_data)["Cluster"] = c("X","X","Y","Y")
collapse_samples(scMethrix_data, colname = "Cluster")
}
