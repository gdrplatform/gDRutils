#' @include global_cache.R

#' Get the identifiers for an identifier field
#'
#' Get the expected identifier(s) for a specific identifier field
#' 
#' @param x string of field to return identifiers for
#' 
#' @return character vector of identifiers for field \code{x}
#' 
#' @details
#' If \code{get_identifier} is called with no arguments, the entire available identifier list is returned.
#' @examples
#' get_identifier("duration")
#' @export
get_identifier <- function(x = NULL) {
  global_cache$get_id(x)
}


#' Set the identifier(s) for an identifier field
#'
#' Set the expected identifier(s) for a specific identifier field
#' 
#' @param x string of field to set identifiers for
#' @param v character vector of identifiers corresponding to identifier field \code{x}
#' 
#' @return a \code{NULL} invisibly 
#' 
#' @details
#' If \code{set_identifier} is called on an identifier that is already set, an error is thrown.
#' @export
set_identifier <- function(x, v) {
  global_cache$set_id(x, v)
  invisible(NULL)
}
