## Identifiers belong to the object. The session list (global_cache.R) is the
## default for new data; a function working on an SE makes that SE's identifiers
## the active ones for its own duration, and restores the session on exit.

#' Use a set of identifiers for the duration of a call
#'
#' \code{local_identifiers()} makes \code{ids} the active identifiers until the
#' calling function returns, on success or on error, then restores the identifiers
#' that were active before. \code{with_identifiers()} does the same around one
#' expression, and \code{with_SE_identifiers()} uses the identifiers stored on a
#' \code{SummarizedExperiment}. Keys absent from \code{ids} keep their current value.
#'
#' Prefer these to \code{set_env_identifier()}, which changes the identifiers for the
#' rest of the session. Called at the top level (the console), there is no function to
#' return from, so \code{local_identifiers()} warns and the identifiers stay set, as with
#' \code{set_env_identifier()}.
#'
#' An SE built by gDRcore stores the identifiers it was built with, and the functions
#' that read an SE (\code{convert_se_assay_to_dt()}, \code{validate_SE()},
#' \code{apply_fit()}, ...) use them in preference to the session's. A key the SE does
#' not store is taken from the session, as \code{get_SE_identifiers()} always did. To
#' correct an object's identifiers, use \code{set_SE_identifiers()}.
#'
#' @param ids named list of identifiers, a subset of \code{names(get_default_identifiers())}.
#' @param .local_envir environment whose exit restores the previous identifiers.
#' @param code expression evaluated with the identifiers active.
#' @param se a \code{SummarizedExperiment}.
#'
#' @return \code{local_identifiers()}: the identifiers active before, invisibly.
#'   \code{with_identifiers()} and \code{with_SE_identifiers()}: the value of \code{code}.
#'
#' @examples
#' f <- function() {
#'   local_identifiers(list(drug = "DrugID"))
#'   get_env_identifiers("drug")
#' }
#' f()                         # "DrugID"
#' get_env_identifiers("drug") # unchanged
#' with_identifiers(list(cellline = "Clid"), get_env_identifiers("cellline"))
#'
#' @keywords identifiers
#' @export
local_identifiers <- function(ids, .local_envir = parent.frame()) {
  checkmate::assert_list(ids, types = "character", names = "unique")
  checkmate::assert_subset(names(ids), names(IDENTIFIERS_LIST))
  checkmate::assert_environment(.local_envir)
  if (identical(.local_envir, globalenv())) {
    warning("local_identifiers() called at the top level: nothing will restore the ",
            "previous identifiers, so this works like set_env_identifier()", call. = FALSE)
  }
  old <- .get_id()
  for (k in names(ids)) {
    .set_id(k, ids[[k]])
  }
  restore <- .restore_ids
  do.call(base::on.exit, list(bquote(.(restore)(.(old))), add = TRUE, after = FALSE),
          envir = .local_envir)
  invisible(old)
}


#' @rdname local_identifiers
#' @export
with_identifiers <- function(ids, code) {
  local_identifiers(ids)
  force(code)
}


#' @rdname local_identifiers
#' @export
with_SE_identifiers <- function(se, code) {
  checkmate::assert_class(se, "SummarizedExperiment")
  local_identifiers(.SE_scope_ids(se))
  force(code)
}


# The identifiers stored on an SE that are known keys with a value. An object without
# stored identifiers, or with only some of them, leaves the rest as they are in the
# session, which is the fallback get_SE_identifiers() has always used. A key that is no
# longer an identifier, or a NULL / non-character value, counts as not stored.
#' @keywords internal
.SE_scope_ids <- function(se) {
  ids <- S4Vectors::metadata(se)[["identifiers"]]
  if (!is.list(ids) || length(ids) == 0L) {
    return(list())
  }
  ids <- ids[intersect(names(ids), names(IDENTIFIERS_LIST))]
  Filter(function(v) is.character(v) && length(v) > 0L, ids)
}


# Entry-point helper: make an SE's identifiers active until the caller returns.
#' @keywords internal
.local_SE_identifiers <- function(se, .local_envir = parent.frame()) {
  local_identifiers(.SE_scope_ids(se), .local_envir = .local_envir)
}


#' @keywords internal
.restore_ids <- function(ids) {
  global_cache$identifiers_list <- ids
  invisible(NULL)
}


#' Canonical form of a set of identifiers, and whether two sets agree
#'
#' \code{canonical_identifiers()} completes \code{ids} with the session's value of
#' every key it lacks (\code{fill}), recording those keys in \code{attr(, "assumed")}:
#' a key an object does not store is read from the session, so it is compared as the
#' session's value too. Keys that are not identifiers (e.g. from an older object) are
#' left out and named in \code{attr(, "unknown")}. Given
#' \code{data}, it also resolves an identifier that lists alternative columns (e.g.
#' \code{barcode = c("Barcode", "Plate")}) to the one present, so two objects built
#' the same way store the same identifiers.
#'
#' \code{same_identifiers()} is \code{TRUE} when two sets name the same columns. A
#' key resolved to one column agrees with a list of alternatives that contains it, and
#' lists of alternatives or tags (\code{barcode}, \code{template},
#' \code{untreated_tag}) are compared as sets; \code{well_position} keeps its order.
#' The keys that differ are in \code{attr(, "differ")}.
#'
#' @param ids,a,b named lists of identifiers.
#' @param data optional \code{data.frame}, or character vector of column names, the
#'   identifiers describe.
#' @param fill named list giving the value of a key \code{ids} lacks. Defaults to the
#'   session's identifiers.
#'
#' @return \code{canonical_identifiers()}: a named list over every identifier key.
#'   \code{same_identifiers()}: a flag.
#'
#' @examples
#' ids <- canonical_identifiers(list(barcode = c("Barcode", "Plate")), data = "Plate")
#' ids$barcode # "Plate"
#' attr(ids, "assumed")
#' same_identifiers(list(barcode = "Plate"), get_default_identifiers())
#'
#' @keywords identifiers
#' @export
canonical_identifiers <- function(ids, data = NULL, fill = get_env_identifiers()) {
  checkmate::assert_list(ids, names = "unique")
  checkmate::assert_list(fill, names = "unique")
  unknown <- setdiff(names(ids), names(IDENTIFIERS_LIST))
  ids <- ids[setdiff(names(ids), unknown)]
  ids <- ids[!vapply(ids, is.null, logical(1))]
  assumed <- setdiff(names(IDENTIFIERS_LIST), names(ids))
  filled <- fill[intersect(assumed, names(fill))]
  out <- c(ids, filled, IDENTIFIERS_LIST[setdiff(assumed, names(fill))])[names(IDENTIFIERS_LIST)]
  if (!is.null(data)) {
    cols <- if (is.character(data)) data else colnames(data)
    checkmate::assert_character(cols)
    for (k in intersect(EXPECT_ONE_IDENTIFIERS, names(out))) {
      if (length(out[[k]]) > 1L) {
        hit <- intersect(out[[k]], cols)
        if (length(hit) == 1L) out[[k]] <- hit
      }
    }
  }
  attr(out, "assumed") <- assumed
  attr(out, "unknown") <- unknown
  out
}


#' @rdname canonical_identifiers
#' @export
same_identifiers <- function(a, b) {
  a <- canonical_identifiers(a)
  b <- canonical_identifiers(b)
  differ <- names(a)[!vapply(names(a), function(k) .same_id_value(k, a[[k]], b[[k]]),
                             logical(1))]
  structure(length(differ) == 0L, differ = differ)
}


# Identifiers whose value is a set (alternative columns, or tags), not an ordered list.
.UNORDERED_IDENTIFIERS <- c("barcode", "template", "untreated_tag")


#' @keywords internal
.same_id_value <- function(k, x, y) {
  if (identical(x, y)) {
    return(TRUE)
  }
  if (k %in% .UNORDERED_IDENTIFIERS && setequal(x, y)) {
    return(TRUE)
  }
  if (k %in% EXPECT_ONE_IDENTIFIERS && xor(length(x) == 1L, length(y) == 1L)) {
    short <- if (length(x) == 1L) x else y
    long <- if (length(x) == 1L) y else x
    return(short %in% long)
  }
  FALSE
}
