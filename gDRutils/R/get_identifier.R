#' @include global_cache.R

#' @export
get_identifier <- function(x = NULL) {
  global_cache$get_id(x)
}


#' @export
set_identifier <- function(x, v) {
  global_cache$set_id(x, v)
}
