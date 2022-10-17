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

#' Get gDR default identifiers required for downstream analysis.
#' @export
get_default_identifiers <- function() {
  IDENTIFIERS_LIST
}

#' Get gDR synonyms for the identifiers
#'
#' Get gDR synonyms for the identifiers
#'
#' @export
get_idfs_synonyms <- function() {
  SYNONYMS_LIST
}

#' Update gDR synonyms for the identifier
#'
#' Update gDR synonyms for the identifier
#'
#' @param data list of charvec with identifiers data
#' @param dict list with dictionary
#'
#' @export
update_idfs_synonyms <- function(data, dict = get_idfs_synonyms()) {

  if (!is.character(data) && !is.list(data)) {
    stop("'data' must be a list of character vector")
  }
  checkmate::assert_list(dict)

  out <- if (is.list(data)) {
    lapply(data, function(x) {
      update_idfs_synonyms(x, dict)
    })
  } else {
    for (idfs in names(dict)) {
      idx <- which(toupper(data) %in% toupper(dict[[idfs]]))
      if (length(idx) > 0) {
        data[idx] <- gDRutils::get_env_identifiers(idfs)
      }
    }
    data
  }
  out
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

#' Update environment identifers from MAE object identifers
#'
#' @param mae_idfs A list containing MAE identifiers
#'
#' @export
update_env_idfs_from_mae <- function(mae_idfs) {
  mae_vs_env <- lapply(mae_idfs, function(x) {
    lapply(names(x), function(i) {
      x[[i]] %in% gDRutils::get_env_identifiers()[[i]]
    })
  })
  for (x in names(mae_idfs)) {names(mae_vs_env[[x]]) <- names(mae_idfs[[x]])}
  if (!all(unlist(mae_vs_env))) {
    # finding the MAE idfs index which needs to be changed
    changed_idfs_exp <- lapply(mae_vs_env, function(x) {
      unlist(unname(lapply(names(x), function(i) {
        if (length(x[[i]]) > 1 && !all(x[[i]])) {
          i
        } else if (length(x[[i]]) == 1 && !all(x[[i]])) {
          i
        }
      })))
    })
    # setting the env identifiers to that of the SE object
    if (length(changed_idfs_exp) > 1L) {
      if (all(changed_idfs_exp[[1]] == changed_idfs_exp[[2]])) {
        lapply(changed_idfs_exp[[1]], function(x) {gDRutils::set_env_identifier(x, mae_idfs[[1]][[x]])})
      } else {
        stop("Both SE objects should have the same identifiers")
      }
    } else {
      lapply(changed_idfs_exp[[1]], function(x) {gDRutils::set_env_identifier(x, mae_idfs[[1]][[x]])})
    }
  }
}
