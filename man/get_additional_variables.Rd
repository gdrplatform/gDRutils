% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_additional_variables}
\alias{get_additional_variables}
\title{Identify and return additional variables in list of dt}
\usage{
get_additional_variables(dt_list, unique = FALSE, prettified = FALSE)
}
\arguments{
\item{dt_list}{list of data.table or data.table containing additional variables}

\item{unique}{logical flag indicating if all variables should be returned
or only those containing more than one unique value}

\item{prettified}{Flag indicating if the provided identifiers in the dt are prettified}
}
\value{
vector of variable names with additional variables
}
\description{
Identify and return additional variables in list of dt
}
\examples{
dt <- data.table::data.table(
  Gnumber = seq_len(10), 
  Concentration = runif(10), 
  Ligand = c(rep(0.5, 5), rep(0, 5))
)
get_additional_variables(dt)

}
\keyword{combination_data}
