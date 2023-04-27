#' @import SummarizedExperiment
#' @importFrom magrittr %$% %>% %<>%
#' @importFrom data.table ":=" .N
#' @importFrom methods is
#' @importFrom tibble tribble

.datatable.aware <- TRUE
NULL

# Prevent R CMD check from complaining about the use of pipe expressions
# standard data.table variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "..req_cols",
      "..opt_cols"
    ),
    utils::packageName())
}
