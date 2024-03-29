% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_duplicated_rows}
\alias{get_duplicated_rows}
\title{Helper function to find duplicated rows}
\usage{
get_duplicated_rows(x, col_names = NULL)
}
\arguments{
\item{x}{data frame}

\item{col_names}{character vector, columns in which duplication are searched for}
}
\value{
integer vector
}
\description{
Helper function to find duplicated rows
}
\examples{
dt <- data.table::data.table(a = c(1, 2, 3), b = c(3, 2, 2))
get_duplicated_rows(dt, "b")
}
\keyword{package_utils}
