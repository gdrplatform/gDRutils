% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duplicates.R
\name{throw_msg_if_duplicates}
\alias{throw_msg_if_duplicates}
\title{throw message if assay data.table contains duplicated rows}
\usage{
throw_msg_if_duplicates(
  dt,
  assay_name = "unknown",
  msg_f = stop,
  preview_max_numb = 4
)
}
\arguments{
\item{dt}{data.table with assay data}

\item{assay_name}{string with the name of the assay}

\item{msg_f}{function to be used to throw the message}

\item{preview_max_numb}{number of rows to preview if duplicates found}
}
\description{
An auxiliary function that checks for duplicated rows in assay data.table,
In case of duplicates it throws a message. The messsage function is by default \code{stop()}
The message function can be customized with \code{msg_f} parameter
}
\examples{
sdata <- get_synthetic_data("finalMAE_small")
smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
throw_msg_if_duplicates(smetrics_data, assay_name = "Metrics", msg_f = futile.logger::flog.info)
}
\keyword{duplicates}
