% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experiment_validators.R
\name{validate_MAE}
\alias{validate_MAE}
\title{Validate MultiAssayExperiment object}
\usage{
validate_MAE(mae)
}
\arguments{
\item{mae}{MultiAssayExperiment object
produced by the gDR pipeline}
}
\value{
\code{NULL} invisibly if the MultiAssayExperiment is valid.
Throws an error if the MultiAssayExperiment is not valid.
}
\description{
Function validates correctness of SE included in MAE by checking multiple cases:
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
validate_MAE(mae)

}
\author{
Bartosz Czech \href{mailto:bartosz.czech@contractors.roche.com}{bartosz.czech@contractors.roche.com}
}
\keyword{experiment}
