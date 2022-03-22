## The following are the default values that will be used for handling MAE experiments
EXPERIMENT_GROUPS <-
  list(`single-agent` = c("single-agent", "cotreatment"),
       matrix = "matrix")

#' get_experiment_groups
#'
#' get experiment groups
#' @param mae MultiAssayExperiment object
#' @export
#'
#' @return list with experiment groups or string (if type not NULL)
#'
#' @author Arkadiusz Gladki <arkadiusz.gladki@@contractors.roche.com>
get_experiment_groups <- function(type = NULL) {
  
  checkmate::assert_choice(type, names(EXPERIMENT_GROUPS), null.ok = TRUE)
  
  if (!is.null(type)) {
    EXPERIMENT_GROUPS[[type]]
  } else {
    EXPERIMENT_GROUPS
  }
  
}
