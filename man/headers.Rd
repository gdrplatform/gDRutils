% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/headers.R
\name{headers}
\alias{headers}
\alias{get_header}
\title{Get or reset headers for one or all header field(s) respectively}
\usage{
get_header(k = NULL)
}
\arguments{
\item{k}{string of field (data type) to return headers for}
}
\value{
For \code{get_header} a character vector of headers for field \code{k}.
}
\description{
Get the expected header(s) for one field or reset all header fields
}
\details{
If \code{get_header} is called with no values, the entire available header list is returned.
}
\examples{
get_header(k = NULL)
get_header("manifest")
}
\keyword{identifiers}
