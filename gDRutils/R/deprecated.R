#' @export
logistic_4parameters <- function(c, Vinf, V0, EC50, h) {
  .Deprecated("predict_efficacy_from_conc")
  predict_efficacy_from_conc(c, Vinf, V0, EC50, h)
}
