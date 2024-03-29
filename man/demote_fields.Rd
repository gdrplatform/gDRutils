% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/concatentate_SEs.R
\name{demote_fields}
\alias{demote_fields}
\title{Demote a metadata field in the \code{rowData} or \code{colData} of a \code{SummarizedExperiment} object
to a nested field of a \code{BumpyMatrix} assay.}
\usage{
demote_fields(se, fields)
}
\arguments{
\item{se}{A \code{SummarizedExperiment} object.}

\item{fields}{Character vector of metadata fields to demote as nested columns.}
}
\value{
A \code{SummarizedExperiment} object with new dimensions resulting from demoting given \code{fields}
to nested columns.
}
\description{
Demote a metadata field in the \code{rowData} or \code{colData} of a \code{SummarizedExperiment} object
to a nested field of a \code{BumpyMatrix} assay.
}
\details{
Revert this operation using \code{promote_fields}.
}
\examples{
mae <- get_synthetic_data("finalMAE_small")
se <- mae[[1]]
se <- promote_fields(se, "ReadoutValue", 2)
demote_fields(se, "ReadoutValue")

}
\seealso{
promote_fields
}
\keyword{SE_operators}
