## The following are the default values that will be used for handling MAE experiments
SUPPORTED_EXPERIMENTS <- c("single-agent", "combination")

#' get_experiment_groups
#'
#' get supported experiments
#' @keywords experiment
#'
#' @return charvec with supported experiment names
#'
#' @examples
#' get_supported_experiments()
#'
#' @export
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
get_supported_experiments <- function() {
  SUPPORTED_EXPERIMENTS
}
