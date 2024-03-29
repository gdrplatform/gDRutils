% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_SE.R
\name{identify_unique_se_metadata_fields}
\alias{identify_unique_se_metadata_fields}
\title{Identify unique metadata fields from a list of \code{SummarizedExperiment}s}
\usage{
identify_unique_se_metadata_fields(SElist)
}
\arguments{
\item{SElist}{named list of \code{SummarizedExperiment}s}
}
\value{
character vector of unique names of metadata
}
\description{
Identify unique metadata fields from a list of \code{SummarizedExperiment}s
}
\examples{
mae <- get_synthetic_data("finalMAE_small")
se <- mae[[1]]
SElist <- list(
  se, 
  se
)
identify_unique_se_metadata_fields(SElist)

}
\keyword{SE_operators}
