## The following are the default values that will be used for handling MAE experiments
EXPERIMENT_GROUPS <-
  list(`single-agent` = c(`single-agent` = "single-agent",
                          `cotreatment` = "cotreatment"),
       matrix = "matrix")

#' get_experiment_groups
#'
#' get experiment groups
#' @param type String indicating the name of an assay group.
#' Defaults to all experiment groups.
#'
#' @return list with experiment groups or string (if type not NULL)
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
