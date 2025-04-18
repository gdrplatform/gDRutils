#' @import SummarizedExperiment
#' @importFrom data.table := .N .SD
#' @importFrom methods is

.datatable.aware <- TRUE

# Prevent R CMD check from complaining the standard data.table variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      ".",
      "GR50",
      "IC50",
      "GR Inf",
      "GR 0",
      "GEC50",
      "h GR",
      "E Inf",
      "E0",
      "EC50",
      "h RV",
      "GR Max",
      "E Max",
      "metric",
      "Response",
      ".N",
      "maxlog10Concentration",
      "MaxEffectiveness",
      "N_conc",
      "normalization_type",
      "RV AOC within set range",
      "GR AOC within set range",
      "project_id",
      "rId",
      "cId",
      "concs",
      "type",
      "name",
      "count",
      "fit_source",
      "summed_conc",
      "LibPath",
      "MaxVersion",
      "Package",
      "UsedVersion",
      "Version"
    ),
    utils::packageName())
}