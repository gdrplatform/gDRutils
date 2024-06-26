#' get_testdata
#'
#' Function to obtain data from gDRtestData and prepare for unit tests
#'
#' @examples 
#' get_testdata()
#'
#' @keywords test_helpers
#' @return list with drugs, cell_lines, raw_data and assay_data
#'
#' @export
get_testdata <- function() {
  
  mae <- get_synthetic_data("finalMAE_small")
  raw_data <- convert_mae_assay_to_dt(mae, "Metrics")
  drug_names <- unique(raw_data$DrugName)
  cell_line_names <- unique(raw_data$CellLineName)
  
  # getting first occurrence of drug_names for each cell_line to avoid aggregation
  dt <- raw_data[, .SD[1], by = c("DrugName", "CellLineName")]
  data.table::setnames(dt,
                       c("CellLineName", "DrugName", "drug_moa", "x_inf", 
                             "x_0", "xc50", "h", "r2", "x_sd_avg", "x_mean", "x_AOC_range", 
                             "x_max", "maxlog10Concentration"),
                       c("Cell Line Name", "Drug Name", "Drug MOA", 
                         "GR Inf", "GR 0", "GEC50", "h GR", "E Inf", "E0", "EC50", "h RV", 
                         "GR Max", "Concentration"))
  dt$GR50 <- dt$EC50
  dt$IC50 <- dt$EC50
  dt$`E Max` <- dt$EC50
  dt$`GR value` <- dt$EC50
  
  list(
    drug_names = drug_names,
    cell_line_names = cell_line_names,
    dt = dt,
    raw_data = raw_data
  )
}

#' get_testdata_combo
#'
#' Function to obtain data from gDRtestData and prepare for unit tests
#'
#' @examples 
#' get_testdata_combo()
#'
#' @keywords test_helpers
#' @return list with drugs, cell_lines, raw_data and assay_data
#'
#' @export
get_testdata_combo <- function() {
  
  mae <- get_synthetic_data("finalMAE_combo_matrix")
  raw_data <- convert_mae_assay_to_dt(mae, "Metrics")
  drug_names <- unique(raw_data$DrugName)
  cell_line_names <- unique(raw_data$CellLineName)
  
  # getting first occurrence of drug_names for each cell_line to avoid aggregation
  dt <- raw_data[, .SD[1], by = c("DrugName", "CellLineName")]
  data.table::setnames(dt,
                       c("CellLineName", "DrugName", "drug_moa", "DrugName_2", "drug_moa_2",
                         "x_inf", "x_0", "xc50", "h", "r2", "x_sd_avg", "x_mean", "x_AOC_range", 
                         "x_max", "maxlog10Concentration"),
                       c("Cell Line Name", "Drug Name", "Drug MOA", "Drug Name 2", "Drug MOA2",
                         "GR Inf", "GR 0", "GEC50", "h GR", "E Inf", "E0", "EC50", "h RV", 
                         "GR Max", "Concentration"))
  dt$GR50 <- dt$EC50
  dt$IC50 <- dt$EC50
  dt$`E Max` <- dt$EC50
  dt$`GR value` <- dt$EC50
  
  list(
    drug_names = drug_names,
    cell_line_names = cell_line_names,
    dt = dt,
    raw_data = raw_data
  )
}

#' get_testdata_codilution
#'
#' Function to obtain data from gDRtestData and prepare for unit tests
#'
#' @examples 
#' get_testdata_codilution()
#'
#' @keywords test_helpers
#' @return list with drugs, cell_lines, raw_data and assay_data
#'
#' @export
get_testdata_codilution <- function() {
  
  mae <- get_synthetic_data("finalMAE_combo_codilution_small")
  raw_data <- convert_mae_assay_to_dt(mae, "Metrics")
  drug_names <- unique(raw_data$DrugName)
  cell_line_names <- unique(raw_data$CellLineName)
  
  # getting first occurrence of drug_names for each cell_line to avoid aggregation
  dt <- raw_data[, .SD[1], by = c("DrugName", "CellLineName")]
  data.table::setnames(dt,
                       c("CellLineName", "DrugName", "drug_moa", "DrugName_2", "drug_moa_2",
                         "x_inf", "x_0", "xc50", "h", "r2", "x_sd_avg", "x_mean", "x_AOC_range", 
                         "x_max", "maxlog10Concentration", "Concentration_2"),
                       c("Cell Line Name", "Drug Name", "Drug MOA", "Drug Name 2", "Drug MOA2",
                         "GR Inf", "GR 0", "GEC50", "h GR", "E Inf", "E0", "EC50", "h RV", 
                         "GR Max", "Concentration", "Concentration 2"))
  dt$GR50 <- dt$EC50
  dt$IC50 <- dt$EC50
  dt$`E Max` <- dt$EC50
  dt$`GR value` <- dt$EC50
  
  list(
    drug_names = drug_names,
    cell_line_names = cell_line_names,
    dt = dt,
    raw_data = raw_data
  )
}

#' gen_synthetic_data
#'
#' Function for generating local synthetic data used for unit tests in modules
#'
#' @param m number of drugs
#' @param n number of records
#' @keywords test_helpers
#'
#' @examples 
#' gen_synthetic_data()
#'
#' @return list with drugs, cell_lines, raw_data and assay_data
#'
#' @export
gen_synthetic_data <- function(m = 1, n = 5) {
  checkmate::assert_number(m, lower = 1)
  checkmate::assert_number(n, lower = 1)
  drug_names <- paste0("drug_00", seq_len(m))
  cell_names <- paste0("cellline_", LETTERS[2:6], "A")
  values <- seq(0.1, 2.5, length.out = m * n)
  dt <- data.table::data.table(
    "Drug Name" = rep(drug_names, n),
    "Drug MOA" = rep(c(rep("moa_A", m - 1), "moa_B"), n),
    "Cell Line Name" = rep(cell_names, each = m),
    "Tissue" = rep(c(rep("tissue_x", n - 1), "tissue_y"), each = m),
    "GR_AOC" = values,
    "GR Inf" = values,
    "GR 0" = values,
    "GEC50" = values,
    "h GR" = values,
    "E Inf" = values,
    "E0" = values,
    "EC50" = values,
    "h RV" = values,
    "GR50" = values,
    "IC50" = values,
    "GR Max" = values,
    "E Max" = values,
    "GR value" = values,
    "Concentration" = values
  )
  list(drug_names = drug_names, cell_names = cell_names, dt = dt)
}
