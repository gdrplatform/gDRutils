% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/identifiers.R
\name{identifiers}
\alias{identifiers}
\alias{get_env_identifiers}
\alias{get_prettified_identifiers}
\alias{set_env_identifier}
\alias{reset_env_identifiers}
\title{Get, set, or reset identifiers for one or all identifier field(s)}
\usage{
get_env_identifiers(k = NULL, simplify = TRUE)

get_prettified_identifiers(k = NULL, simplify = TRUE)

set_env_identifier(k, v)

reset_env_identifiers()
}
\arguments{
\item{k}{String corresponding to identifier name.}

\item{simplify}{Boolean indicating whether output should be simplified.}

\item{v}{Character vector corresponding to the value for given identifier \code{k}.}
}
\value{
For any \code{set}ting or \code{reset}ting functionality, a \code{NULL} invisibly.
For \code{get_env_identifiers} a character vector of identifiers for field \code{k}.
For functions called with no arguments, the entire available identifier list is returned.

list or charvec depends on unify param

list or charvec depends on unify param

\code{NULL}

\code{NULL}
}
\description{
Get, set, or reset the expected identifier(s) for one or all identifier field(s).
Identifiers are used by the gDR processing functions to identify which columns in a \code{data.table}
correspond to certain expected fields. Functions of the family \code{*et_identifier} will look for
identifiers from the environment while functions of the family \code{*et_SE_identifiers} will look for
identifiers in the \code{metadata} slot of a \code{SummarizedExperiment} object.
See details for expected identifiers and their definitions.
}
\details{
Identifiers supported by the gDR suite include:
\itemize{
\item{"barcode": String of column name containing barcode metadata}
\item{"cellline": String of column name containing unique, machine-readable cell line identifiers}
\item{"cellline_name": String of column name containing human-friendly cell line names}
\item{"cellline_tissue": String of column name containing metadata on cell line tissue type}
\item{"cellline_ref_div_time": String of column name containing reference division time for cell lines}
\item{"cellline_parental_identifier": String of column name containing unique, machine-readable
parental cell line identifiers. Used in the case of derived or engineered cell lines.}
\item{"drug": String of column name containing unique, machine-readable drug identifiers}
\item{"drug_name": String of column name containing human-friendly drug names}
\item{"drug_moa": String of column name containing metadata for drug mode of action}
\item{"duration": String of column name containing metadata on duration that cells were treated (in hours)}
\item{"template": String of column name containing template metadata}
\item{"untreated_tag": Character vector of entries that identify control, untreated wells}
\item{"well_position": Character vector of column names containing metadata on well positions on a plate}
}
}
\examples{
get_env_identifiers("duration") # "Duration"

}
\keyword{identifiers}
