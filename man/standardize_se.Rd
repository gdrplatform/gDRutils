% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize_MAE.R
\name{standardize_se}
\alias{standardize_se}
\title{Standardize SE by switching from custom identifiers into gDR-default}
\usage{
standardize_se(se, use_default = TRUE)
}
\arguments{
\item{se}{a SummarizedExperiment object with drug-response data generate by gDR pipeline}

\item{use_default}{boolean indicating whether or not to use default
identifiers for standardization}
}
\value{
se a SummarizedExperiment with default gDR identifiers
}
\description{
Standardize SE by switching from custom identifiers into gDR-default
}
\examples{
mae <- get_synthetic_data("finalMAE_small.qs")
se <- mae[[1]]
S4Vectors::metadata(se)$identifiers$drug <- "druug"
standardize_se(se)

}
\keyword{standardize_MAE}
