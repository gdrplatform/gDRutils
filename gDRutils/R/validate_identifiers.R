#' Check that specified identifier values exist in the data.
#' 
#'
#'
#' @param df data.frame with \code{colnames}.
#' @param identifiers Named list of identifiers.
#' If not passed, defaults to \code{get_env_identifiers()}.
#'
#' @return Named list of identifiers.
#'
#' @details
#' Note that this does NOT set the identifiers anywhere (i.e. environment or \code{SummarizedExperiment} object).
#' If identifiers do not validate, will throw error as side effect.
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

  identifiers <- .singlefy_polymapped_identifiers(df, exp_one_ids, id_map)
  .check_identifiers(df, identifiers, exp_one_ids, req_ids)

  identifiers
}


#' Ensure all identifiers that are expected to have only one identifier indeed have only one.
#'
#' In some cases, the environment identifiers may have a many-to-one mapping between
#' the standardized identifier and its values. In such cases, this function will check
#' that only the identifier that is present in the data is in the returned identifiers.
#'
#' @param df data.frame with named columns
#' @param exp_one_ids Character vector of identifiers where only one column name is expected. 
#' @param id_map Named list where values represent column names in \code{df} and
#' \code{names()} of list represent standardized identifier names.
#'
#' @return String describing validation failures. \code{NULL} if no failures.
#' 
#' @details
#' Note that this will not change the environment identifiers.
#' @noRd
.singlefy_polymapped_identifiers <- function(df, exp_one_ids, id_map) {
  exp_one_id_maps <- id_map[names(id_map) %in% exp_one_ids]
  polymappings <- exp_one_id_maps[lengths(exp_one_id_maps) > 1L]

  npoly <- length(polymappings)
  if (npoly > 0L) {
    for (id in names(polymappings)) {
      shared <- intersect(polymappings[[id]], colnames(df))
      n_shared <- length(shared)
      new_v <- if (nshared == 0L) {
        NA
      } else if (nshared == 1L) {
        shared
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
#' @param req_ids Character vector of identifiers required in the \code{data.frame} 
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


#' @return String or \code{NULL} of message
#' @noRd
.check_polymapped_identifiers <- function(df, exp_one_ids, id_map) {
  exp_one_id_maps <- id_map[names(id_map) %in% exp_one_ids]
  polymappings <- exp_one_id_maps[lengths(exp_one_id_maps) > 1L]

  msg <- NULL
  if (length(polymappings) > 0L) {
    multi <- vapply(polymappings, function(x) {length(intersect(polymappings[[x]], colnames(df))) > 1}, boolean(0))
    msg <- paste0("more than one mapping for identifier(s):", paste0(names(polymappings)[multi], collapse = ", "), "\n")
  }
  msg
}
