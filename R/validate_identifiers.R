#' Check that specified identifier values exist in the data.
#' 
#' Check that specified identifier values exist in the data and error otherwise.
#'
#' @param df data.frame with \code{colnames}.
#' @param identifiers Named list of identifiers where the \code{names} are standardized identifier names.
#' If not passed, defaults to \code{get_env_identifiers()}.
#' @param req_ids Character vector of standardized identifier names required to pass identifier validation.
#' @param exp_one_ids Character vector of standardized identifiers names
#' where only one identifier value is expected. 
#' If not passed, defaults to \code{get_expect_one_identifiers()}.
#'
#' @return Named list of identifiers modified to pass validation against the input data.
#' Errors with explanatory message if validation cannot pass with the given identifiers and data.
#'
#' @details
#' Note that this does NOT set the identifiers anywhere (i.e. environment or \code{SummarizedExperiment} object).
#' If identifiers do not validate, will throw error as side effect.
#' 
#' @examples 
#' validate_identifiers(
#'   DataFrame("Barcode" = NA, "Duration" = NA, "Template" = NA, "clid" = NA), 
#'   req_ids = "barcode"
#' )
#' 
#' @export
validate_identifiers <- function(df, identifiers = NULL, req_ids = NULL, exp_one_ids = NULL) {
  if (is.null(identifiers)) {
    identifiers <- get_env_identifiers()
  }
  if (is.null(req_ids)) {
    req_ids <- get_required_identifiers()
  }
  if (is.null(exp_one_ids)) {
    exp_one_ids <- get_expect_one_identifiers()
  }

  identifiers <- .modify_polymapped_identifiers(df, exp_one_ids, identifiers)
  .check_identifiers(df, identifiers, exp_one_ids, req_ids)

  identifiers
}


#' Modify identifier values to reflect the data.
#'
#' Identifier mappings may have a many-to-one relationship between
#' standardized identifiers and their values. This many-to-one relationship is used to store
#' "default values" in cases where a column may be named in many ways depending on the source. 
#' In such cases, this function will replace the
#' many-to-one mapping with only the identifier value that is present in the data.
#'
#' @param df data.frame with named columns
#' @param exp_one_ids Character vector of standardized identifiers
#' where only one identifier value is expected. 
#' @param id_map Named list where values represent column names in \code{df} and
#' \code{names()} of list represent standardized identifier names.
#'
#' @return Named list of modified identifiers for many-to-one mappings.
#' 
#' @details
#' This is most often used when the identifiers are set to include some default values.
#' Note that this will not change the environment identifiers.
#' @noRd
.modify_polymapped_identifiers <- function(df, exp_one_ids, id_map) {
  exp_one_id_maps <- id_map[names(id_map) %in% exp_one_ids]
  polymappings <- exp_one_id_maps[lengths(exp_one_id_maps) > 1L]

  npoly <- length(polymappings)
  if (npoly > 0L) {
    for (id in names(polymappings)) {
      shared <- intersect(polymappings[[id]], colnames(df))
      nshared <- length(shared)
      if (nshared == 0L) {
        new_v <- NA
      } else if (nshared == 1L) {
        new_v <- shared
      }  else {
        stop(sprintf("multiple valid identifier values found for identifier: '%s'", id))
      }
      id_map[[id]] <- new_v
    }
  }

  return(id_map)
}


#' @noRd
.check_identifiers <- function(df, id_map, exp_one_ids, req_ids) {
  msg <- NULL
  msg1 <- .check_polymapped_identifiers(df, exp_one_ids, id_map)
  msg2 <- .check_required_identifiers(df, req_ids, id_map)
  msg <- c(msg, msg1, msg2)
  if (!is.null(msg)) {
    stop(msg)
  }
  invisible(NULL)
}


#' Ensure all required identifiers are present in the data.
#'
#' Ensure all required identifiers are present in the data according to an identifier map.
#'
#' @param df data.frame with named columns
#' @param req_ids Character vector of identifiers required to pass identifier validation.
#' @param id_map Named list where values represent column names in \code{df} and
#' \code{names()} of list represent standardized identifier names.
#'
#' @return String describing validation failures. \code{NULL} if no failures.
#' @noRd
.check_required_identifiers <- function(df, req_ids, id_map) {
  missing <- !req_ids %in% names(id_map)
  if (any(missing)) {
    stop(sprintf("required identifiers: '%s' missing in 'id_map'",
      paste0(req_ids[missing], collapse = ", ")))
  }

  gt_one <- lengths(id_map[req_ids]) != 1L
  if (any(gt_one)) {
    stop(sprintf("more than one identifier value found for required identifiers: '%s'",
      paste0(names(id_map[req_ids][gt_one]), collapse = ", ")))
  }

  msg <- NULL
  req_map <- id_map[req_ids]
  present <- unlist(unname(req_map)) %in% colnames(df)
  if (!all(present)) {
    err_map <- req_map[!present]
    msg <- sprintf("specified value identifier(s): '%s' do not exist for standardized identifier(s): '%s'\n",
      paste(unname(err_map), sep = ", "),
      paste(names(err_map), sep = ", "))
  }
  msg
}


#' @return String of message or \code{NULL}.
#' @noRd
.check_polymapped_identifiers <- function(df, exp_one_ids, id_map) {
  exp_one_id_maps <- id_map[names(id_map) %in% exp_one_ids]
  polymappings <- exp_one_id_maps[lengths(exp_one_id_maps) > 1L]

  msg <- NULL
  if (length(polymappings) > 0L) {
    msg <- sprintf("more than one mapping for identifier(s): '%s'\n", paste0(names(polymappings), collapse = ", "))
  }
  msg
}
