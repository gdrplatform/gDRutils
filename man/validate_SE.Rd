% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experiment_validators.R
\name{validate_SE}
\alias{validate_SE}
\title{Validate SummarizedExperiment object}
\usage{
validate_SE(se, expect_single_agent = FALSE)
}
\arguments{
\item{se}{SummarizedExperiment object
produced by the gDR pipeline}

\item{expect_single_agent}{a logical indicating if the function
should check whether the SummarizedExperiment is single-agent data}
}
\value{
\code{NULL} invisibly if the SummarizedExperiment is valid.
Throws an error if the SummarizedExperiment is not valid.
}
\description{
Function validates correctness of SE by checking multiple cases:
\itemize{
\item detection of duplicated rowData/colData,
\item incompatibility of rownames/colnames,
\item occurrence of necessary assays,
\item detection of mismatch of CLIDs inside colData and colnames (different order),
\item correctness of metadata names.
}
}
\examples{
mae <- get_synthetic_data("finalMAE_small")
se <- mae[[1]]
validate_SE(se)

}
\keyword{experiment}
