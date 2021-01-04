#' @include global_cache.R

#' @title Get, set, or reset identifiers for one or all identifier field(s)
#'
#' @description Get, set, or reset the expected identifier(s) for one or all identifier field(s)
#'
#' @return 
#' For \code{set_identifier} or \code{reset_identifiers} a \code{NULL} invisibly.
#' For \code{get_identifier} a character vector of identifiers for field \code{k}.
#' @examples
#' \dontrun{
#' get_identifier("duration")
#' set_identifier("duration", "SET DURATION TIME")
#' get_identifier("duration")
#' reset_identifiers() 
#' get_identifier("duration")
#'}
#' @name identifiers
NULL


#' @param k string of field to get or set identifiers for
#' @details
#' If \code{get_identifier} is called with no arguments, the entire available identifier list is returned.
#' @rdname identifiers
#' @export
#' 
get_identifier <- function(k = NULL) {
  global_cache$get_id(k)
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
  global_cache$set_id(k, v)
  invisible(NULL)
}


#' @details
#' \code{reset_identifiers} should be used if processing multiple files within the same R session
#' that use different file identifiers.  
#' @rdname identifiers
#' @export
#' 
reset_identifiers <- function() {
  global_cache$reset_ids()
  invisible(NULL)
}
