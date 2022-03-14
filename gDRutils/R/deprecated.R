#' @export
logistic_4parameters <- function(c, Vinf, V0, EC50, h) {
  .Deprecated("predict_efficacy_from_conc")
  predict_efficacy_from_conc(c, Vinf, V0, EC50, h)
}
  

#' @export
get_identifier <- function(k = NULL, simplify = TRUE) {
  .Deprecated("get_env_identifiers")
  get_env_identifiers(k = k, simplify = simplify)
}


#' @export
set_identifier <- function(k, v) {
  .Deprecated("set_env_identifier")
  set_env_identifier(k = k, v = v)
}


#' @export
reset_identifier <- function(k, v) {
  .Deprecated("reset_env_identifiers")
  reset_env_identifiers()
}
