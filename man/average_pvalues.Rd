% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{average_pvalues}
\alias{average_pvalues}
\title{Average p-values using Fisher's method
Combines a vector of p-values into a single representative p-value.
It implements Fisher's method, where the test statistic is calculated as
\deqn{X_{2k}^2 = -2 \sum_{i=1}^{k} \ln(p_i)}.
This statistic follows a chi-squared distribution with 2k degrees of freedom (where k is the number
of p-values), from which the combined p-value is derived.}
\usage{
average_pvalues(p_values)
}
\arguments{
\item{p_values}{A numeric vector of p-values. Values are expected to be between 0 and 1.
The function assumes at least one non-NA value is provided.}
}
\value{
A single, combined p-value as a numeric value.
}
\description{
Average p-values using Fisher's method
Combines a vector of p-values into a single representative p-value.
It implements Fisher's method, where the test statistic is calculated as
\deqn{X_{2k}^2 = -2 \sum_{i=1}^{k} \ln(p_i)}.
This statistic follows a chi-squared distribution with 2k degrees of freedom (where k is the number
of p-values), from which the combined p-value is derived.
}
\keyword{internal}
