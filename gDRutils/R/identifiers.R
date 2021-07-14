#' @title Get, set, or reset identifiers for one or all identifier field(s)
#'
#' @description Get, set, or reset the expected identifier(s) for one or all identifier field(s).
#' Identifiers are used by the gDR processing functions to identify which columns in a \code{data.frame}
#' correspond to certain expected fields. See details for expected identifiers.
#'
#' @param k String corresponding to identifier name.
#' @param v Character vector corresponding to the value for given identifier \code{k}.
#'
#' @return 
#' For \code{set_identifier} or \code{reset_identifiers} a \code{NULL} invisibly.
#' For \code{get_identifier} a character vector of identifiers for field \code{k}.
#'
#' @examples
#' \dontrun{
#' get_identifier("duration") # "Duration"
#' set_identifier("duration", "Duration_Time")
#' get_identifier("duration") # "Duration_Time"
#' reset_identifiers() 
#' get_identifier("duration") # "Duration"
#'}
#'
#' @details
#' Identifiers supported by the gDR suite include:
#' \itemize{
#'  \item{"barcode": }{String of column name containing barcode metadata}
#'  \item{"cellline": }{String of column name containing unique, machine-readable cell line identifiers}
#'  \item{"cellline_name": }{String of column name containing human-friendly cell line names}
#'  \item{"cellline_tissue": }{String of column name containing metadata on cell line tissue type}
#'  \item{"cellline_ref_div_time": }{String of column name containing reference division time for cell lines}
#'  \item{"cellline_parental_identifier": }{String of column name containing unique, machine-readable 
#'    parental cell line identifiers. Used in the case of derived or engineered cell lines.}
#'  \item{"drug": }{String of column name containing unique, machine-readable drug identifiers}
#'  \item{"drugname": }{String of column name containing human-friendly drug names}
#'  \item{"drug_moa": }{String of column name containing metadata for drug mode of action}
#'  \item{"duration": }{String of column name containing metadata on duration that cells were treated (in hours)}
#'  \item{"template": }{String of collumn name containing template metadata}
#'  \item{"untreated_tag": }{Character vector of entries that identify control, untreated wells}
#'  \item{"well_position": }{Character vector of column names containing metadata on well positions on a plate}
#' }
#'
#' @name identifiers
NULL


#' @param k string of field to get or set identifiers for
#' @details
#' If \code{get_identifier} is called with no arguments, the entire available identifier list is returned.
#' @rdname identifiers
#' @export
#' 
get_identifier <- function(k = NULL) {
  .get_id(k)
}


#' @param k string of field to get or set identifiers for
#' @param v character vector of identifiers corresponding to identifier field \code{k}
#' @details
#' \code{set_identifier} can be called on an identifier that is already set.
#' A common use case to do so is when trying to load multiple files of different identifiers
#' in the same R session.
#' @rdname identifiers
#' @export
#' 
set_identifier <- function(k, v) {
  .set_id(k, v)
}


#' @details
#' \code{reset_identifiers} should be used if processing multiple files within the same R session
#' that use different file identifiers.  
#' @rdname identifiers
#' @export
#' 
reset_identifiers <- function() {
  .reset_ids()
}
