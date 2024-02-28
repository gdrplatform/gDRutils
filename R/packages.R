#' @import SummarizedExperiment
#' @importFrom data.table := .N .SD
#' @importFrom methods is

.datatable.aware <- TRUE



# Prevent R CMD check from complaining about the use of pipe expressions
# standard data.table variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "rId",
      "cId",
      "concs",
      "type",
      "name"
    ), 
    utils::packageName())
}
