#' @title Get, set, or reset identifiers for one or all identifier field(s)
#'
#' @description Get, set, or reset the expected identifier(s) for one or all identifier field(s).
#' Identifiers are used by the gDR processing functions to identify which columns in a \code{data.frame}
#' correspond to certain expected fields. Functions of the family \code{*et_identifier} will look for 
#' identifiers from the environment while functions of the family \code{*et_SE_identifiers} will look for
#' identifiers in the \code{metadata} slot of a \code{SummarizedExperiment} object.
#' See details for expected identifiers and their definitions.
#'
#' @param k String corresponding to identifier name.
#' @param v Character vector corresponding to the value for given identifier \code{k}.
#' @param simplify Boolean indicating whether output should be simplified.
#'
#' @return 
#' For any \code{set}ting or \code{reset}ting functionality, a \code{NULL} invisibly.
#' For \code{get_env_identifiers} a character vector of identifiers for field \code{k}.
#' For functions called with no arguments, the entire available identifier list is returned.
#'
#' @examples
#' \dontrun{
#' get_env_identifiers("duration") # "Duration"
#' set_env_identifier("duration", "Duration_Time")
#' get_env_identifiers("duration") # "Duration_Time"
#' reset_env_identifiers() 
#' get_env_identifiers("duration") # "Duration"
#'}
#'
#' @details
#' Identifiers supported by the gDR suite include:
#' \itemize{
#'  \item{"barcode": String of column name containing barcode metadata}
#'  \item{"cellline": String of column name containing unique, machine-readable cell line identifiers}
#'  \item{"cellline_name": String of column name containing human-friendly cell line names}
#'  \item{"cellline_tissue": String of column name containing metadata on cell line tissue type}
#'  \item{"cellline_ref_div_time": String of column name containing reference division time for cell lines}
#'  \item{"cellline_parental_identifier": String of column name containing unique, machine-readable 
#'    parental cell line identifiers. Used in the case of derived or engineered cell lines.}
#'  \item{"drug": String of column name containing unique, machine-readable drug identifiers}
#'  \item{"drug_name": String of column name containing human-friendly drug names}
#'  \item{"drug_moa": String of column name containing metadata for drug mode of action}
#'  \item{"duration": String of column name containing metadata on duration that cells were treated (in hours)}
#'  \item{"template": String of collumn name containing template metadata}
#'  \item{"untreated_tag": Character vector of entries that identify control, untreated wells}
#'  \item{"well_position": Character vector of column names containing metadata on well positions on a plate}
#' }
#'
#' @name identifiers
NULL


#' @rdname identifiers
#' @export
#' 
get_env_identifiers <- function(k = NULL, simplify = TRUE) {
  if (simplify) {
    if (length(k) > 1L) {
      stop("more than one identifier found, please use argument: `simplify = FALSE`")
    } else {
      .get_id(k)
    }
  } else {
    id_vector <- Vectorize(function(i) .get_id(i), SIMPLIFY = FALSE)
    id_vector(k)
  }
}


#' @rdname identifiers
#' @export
get_prettified_identifiers <- function(k = NULL, simplify = TRUE) {
  idfs <- get_env_identifiers(k, simplify = simplify)
  pidfs <- prettify_flat_metrics(idfs, human_readable = TRUE)
  if (is.null(k)) {
    names(pidfs) <- names(idfs)
    # dirty hack for 'untreated_tag' and 'well_position' which are charvec improperly prettified
    # TODO: fix this issue on the prettify_* function level
    pidfs["untreated_tag"] <- idfs["untreated_tag"]
    pidfs["well_position"] <- idfs["well_position"]
  }
  pidfs
}


#' Get identifiers required for downstream analysis.
#' @export
get_required_identifiers <- function() {
  REQ_COL_IDENTIFIERS
}


#' Get identifiers that expect only one value for each identifier.
#' @export
get_expect_one_identifiers <- function() {
  EXPECT_ONE_IDENTIFIERS
}


#' @rdname identifiers
#' @export
#' 
set_env_identifier <- function(k, v) {
  .set_id(k, v)
}


#' @rdname identifiers
#' @export
#' 
reset_env_identifiers <- function() {
  .reset_ids()
}

#' Get descriptions for identifiers
#'
#' @param k identifier key
#'
#' @export
get_identifers_desc <- function(k = NULL) {
  checkmate::assert_string(k, null.ok = TRUE)
  desc <- yaml::read_yaml(system.file(package = "gDRutils", "identifier_descriptions.yaml"))
  if (is.null(k)) {
    desc
  } else {
    checkmate::assert_true(k %in% names(desc))
    desc[[k]]
  }
}
