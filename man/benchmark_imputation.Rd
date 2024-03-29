% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scMethrix_plot.R
\name{benchmark_imputation}
\alias{benchmark_imputation}
\title{Evaluates imputations methods by NRMSE or AUC}
\usage{
benchmark_imputation(
  scm = NULL,
  assay = "score",
  sparse_prop = seq(0.1, 0.9, 0.1),
  iterations = 3,
  imp_methods = c(iPCA = function(...) impute_regions(type = "iPCA", ...), RF =
    function(...) impute_regions(type = "RF", ...), kNN = function(...)
    impute_regions(type = "kNN", ...)),
  type = "RMSE"
)
}
\arguments{
\item{scm}{\code{\link{scMethrix}}; the single cell methylation experiment}

\item{assay}{string; name of an existing assay. Default = "score"}

\item{sparse_prop}{numeric; A sparsity proportion between 0 and 1. E.g. 0.1 replaces 10\% of the matrix with NA}

\item{iterations}{integer; Number of iterations to test}

\item{imp_methods}{closure; The imputation methods to compare.}

\item{type}{character; descriptive statistic. Can be either "AUC" or "RMSE". Default "RMSE"}
}
\value{
ggplot; The graph showing the NRMSE for each imputation method at each sparsity
}
\description{
Evaluates imputations methods by NRMSE or AUC
}
\details{
Does stuff
}
\examples{
data('scMethrix_data')
\dontrun{
scMethrix_data <- impute_regions(scMethrix_data, new_assay="impute",type="RF")
benchmark_imputation(scMethrix_data, assay="impute", sparse_prop = c(0.1,0.5,0.85))
}
}
