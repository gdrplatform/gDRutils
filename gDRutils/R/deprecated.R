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
