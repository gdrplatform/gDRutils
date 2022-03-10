#' Check that specified identifiers exist in the data.
#' 
#' @param df data.frame with \code{colnames}.
#' @param identifiers Named list of identifiers.
#' If not passed, defaults to \code{get_env_identifiers}.
#'
#' @return Named list of identifiers.
#'
#' @details
#' Note that this does NOT set the identifiers.
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
    req_ids <- get_expect_one_identifiers()
  }

  msg <- NULL
  msg1 <- singlefy_polymapped_identifiers(df, exp_one_ids, id_map)
  # TODO: This is too complicated. Fix me. 
  if (is.list(msg1)) {
    identifiers <- msg1
    msg1 <- NULL
  }
  msg2 <- check_required_identifiers(df, req_ids, identifiers)

  msg <- c(msg, msg1, msg2)
  if (!is.null(msg)) {
    stop(msg)
  }
  identifiers
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
#' @export
check_required_identifiers <- function(df, req_ids, id_map) {
  # TODO: check only one value for each identifier.
  msg <- NULL
  req_map <- id_map[req_ids]
  present <- unlist(unname(req_map)) %in% colnames(df)
  if (!all(present)) {
    err_map <- req_map[!present]
    msg <- sprintf("specified value identifier(s): '%s' do not exist for standardized identifier(s): '%s'",
      paste(unname(err_map), sep = ", "),
      paste(names(err_map), sep = ", "))
  }
  msg
}

#' Ensure all 
#' @param df data.frame with named columns
#' @param exp_one_ids Character vector of identifiers where only one column name is expected. 
#' @param id_map Named list where values represent column names in \code{df} and
#' \code{names()} of list represent standardized identifier names.
#'
#' @return String describing validation failures. \code{NULL} if no failures.
#' @export
singlefy_polymapped_identifiers <- function(df, exp_one_ids, id_map) {
  exp_one_id_maps <- id_map[names(id_map) %in% exp_one_ids]
  polymappings <- exp_one_id_maps[lengths(exp_one_id_maps) > 1L]

  npoly <- length(polymappings)
  msg <- NULL
  if (npoly > 0L) {
    for (id in names(polymappings)) {
      shared <- intersect(polymappings[[id]], colnames(df))
      n_shared <- length(shared)
      if (nshared > 1L) {
        msg <- c(msg, paste0("more than one mapping for identifier:", id, "\n"))
      } else {
        new_v <- if (nshared == 0L) {
          NA
        } else if (nshared == 1L) {
          shared
        } 
        id_map[[id]] <- new_v
      }
    }
  }

  # TODO: Fxn should have only one job. Fix me.
  if (is.null(msg)) {
    return(id_map)
  } else {
    return(msg)
  }
}
