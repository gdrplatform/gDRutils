#' Deprecated function, use predict_efficacy_from_conc.
#' @param c Numeric vector representing concentrations to predict efficacies for.
#' @param Vinf Numeric vector representing the asymptotic value of the sigmoidal fit to the dose-response
#'  data as concentration goes to infinity.
#' @param V0 Numeric vector representing the asymptotic metric value corresponding to a concentration of 0
#'  for the primary drug.
#' @param EC50 Numeric vector representing the drug concentration at half-maximal effect.
#' @param h Numeric vector representing the hill coefficient of the fitted curve, which reflects how steep
#'  the dose-response curve is.
#' @export
logistic_4parameters <- function(c, Vinf, V0, EC50, h) {
  .Deprecated("predict_efficacy_from_conc")
  predict_efficacy_from_conc(c, Vinf, V0, EC50, h)
}
  
#' Deprecated function, use get_env_identifier.
#' @param k String corresponding to identifier name.
#' @param simplify Boolean indicating whether output should be simplified.
#' @export
get_identifier <- function(k = NULL, simplify = TRUE) {
  .Deprecated("get_env_identifiers")
  get_env_identifiers(k = k, simplify = simplify)
}


#' Deprecated function, use set_env_identifier.
#' @param k String corresponding to identifier name.
#' @param v Character vector corresponding to the value for given identifier \code{k}.
#' @export
set_identifier <- function(k, v) {
  .Deprecated("set_env_identifier")
  set_env_identifier(k = k, v = v)
}


#' Deprecated function, use reset_env_identifier.
#' @param k String corresponding to identifier name.
#' @param v Character vector corresponding to the value for given identifier \code{k}.
#' @export
reset_identifier <- function(k, v) {
  .Deprecated("reset_env_identifiers")
  reset_env_identifiers()
}
