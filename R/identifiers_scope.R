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
#' rest of the session.
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
  checkmate::assert_list(ids, names = "unique")
  checkmate::assert_subset(names(ids), names(IDENTIFIERS_LIST))
  checkmate::assert_environment(.local_envir)
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


# The identifiers stored on an SE that are known keys. An object without stored
# identifiers, or with only some of them, leaves the rest as they are in the session,
# which is the fallback get_SE_identifiers() has always used.
#' @keywords internal
.SE_scope_ids <- function(se) {
  ids <- S4Vectors::metadata(se)[["identifiers"]]
  if (!is.list(ids) || length(ids) == 0L) {
    return(list())
  }
  ids[intersect(names(ids), names(IDENTIFIERS_LIST))]
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
#' \code{canonical_identifiers()} completes \code{ids} with the default value of
#' every key it lacks, recording those keys in \code{attr(, "assumed")}. Given
#' \code{data}, it also resolves an identifier that lists alternative columns (e.g.
#' \code{barcode = c("Barcode", "Plate")}) to the one present, so two objects built
#' the same way store the same identifiers.
#'
#' \code{same_identifiers()} is \code{TRUE} when two sets name the same columns. A
#' key resolved to one column agrees with a list of alternatives that contains it.
#' The keys that differ are in \code{attr(, "differ")}.
#'
#' @param ids,a,b named lists of identifiers.
#' @param data optional \code{data.frame}, or character vector of column names, the
#'   identifiers describe.
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
canonical_identifiers <- function(ids, data = NULL) {
  checkmate::assert_list(ids, names = "unique")
  unknown <- setdiff(names(ids), names(IDENTIFIERS_LIST))
  if (length(unknown) > 0L) {
    stop(sprintf("unknown identifier key(s): %s", toString(sQuote(unknown, FALSE))))
  }
  assumed <- setdiff(names(IDENTIFIERS_LIST), names(ids))
  out <- c(ids, IDENTIFIERS_LIST[assumed])[names(IDENTIFIERS_LIST)]
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


#' @keywords internal
.same_id_value <- function(k, x, y) {
  if (identical(x, y)) {
    return(TRUE)
  }
  if (k %in% EXPECT_ONE_IDENTIFIERS && xor(length(x) == 1L, length(y) == 1L)) {
    short <- if (length(x) == 1L) x else y
    long <- if (length(x) == 1L) y else x
    return(short %in% long)
  }
  FALSE
}
