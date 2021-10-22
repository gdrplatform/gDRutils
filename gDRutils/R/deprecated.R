#' @export
logistic_4parameters <- function(c, Vinf, V0, EC50, h) {
  .Deprecated("evaluate_efficacy_from_conc")
  evaluate_efficacy_from_conc(c, Vinf, V0, EC50, h)
}
