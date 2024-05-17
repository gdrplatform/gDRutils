## The following are the default values that will be used for handling MAE experiments
SUPPORTED_EXPERIMENTS <-
  c(sa = "single-agent", combo = "combination")
EXPERIMENT_GROUPS <-
  list(`single-agent` = c(`single-agent` = "single-agent",
                          `co-dilution` = "co-dilution"),
       combination = "combination")

#' get_supported_experiments
#'
#' get supported experiments
#' @keywords experiment
#
#' @param type String indicating the type of experiment
#
#' @return charvec with supported experiment name(s)
#'
#' @examples
#' get_supported_experiments()
#'
#' @export
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
get_supported_experiments <- function(type = NULL) {
  
  checkmate::assert_choice(type, names(SUPPORTED_EXPERIMENTS), null.ok = TRUE)
  
  if (is.null(type)) {
    as.character(SUPPORTED_EXPERIMENTS)
  } else {
    as.character(SUPPORTED_EXPERIMENTS[type])
  }
}

#' get_experiment_groups
#'
#' get experiment groups
#' @param type String indicating the name of an assay group.
#' Defaults to all experiment groups.
#' @keywords experiment
#'
#' @return list with experiment groups or string (if type not NULL)
#'
#' @examples 
#' get_experiment_groups()
#' 
#' @export
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
get_experiment_groups <- function(type = NULL) {
  
  checkmate::assert_choice(type, names(EXPERIMENT_GROUPS), null.ok = TRUE)
  
  if (!is.null(type)) {
    EXPERIMENT_GROUPS[[type]]
  } else {
    EXPERIMENT_GROUPS
  }
  
}
