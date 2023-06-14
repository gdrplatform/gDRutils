#' gen_testdata
#'
#' Function to obtain data from gDRtestData and prepare for unit tests
#'
#' @examples 
#' gen_testdata()
#'
#' @return list with drugs, cell_lines, raw_data and assay_data
#'
#' @export
gen_testdata <- function() {
  
  mae <- readRDS(system.file("testdata", "finalMAE_small.RDS", package = "gDRtestData"))
  raw_data <- gDRutils::convert_mae_assay_to_dt(mae, "Metrics")
  drug_names <- unique(raw_data$DrugName)
  cell_line_names <- unique(raw_data$CellLineName)
  
  # getting first occurrence of drug_names for each cell_line to avoid aggregation
  increment <- table(raw_data$DrugName)[1]
  idx <- rep(match(cell_line_names, raw_data$CellLineName), times = length(drug_names)) +
    rep(seq(0, nrow(raw_data) - increment, increment), each = length(drug_names))
  
  dt <- raw_data[idx, ]
  data.table::setnames(dt,
                       c("CellLineName", "DrugName", "Tissue", "drug_moa", "x_inf", 
                             "x_0", "xc50", "h", "r2", "x_sd_avg", "x_mean", "x_AOC_range", 
                             "x_max", "maxlog10Concentration"),
                       c("Cell Line Name", "Drug Name", "Primary Tissue", "Drug MOA", 
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




#' gen_synthetic_data
#'
#' Function for generating local synthetic data used for unit tests in modules
#'
#' @examples 
#' gen_synthetic_data()
#'
#' @return list with drugs, cell_lines, raw_data and assay_data
#'
#' @export
gen_synthetic_data <- function(m, n) {
  drug_names <- paste0("drug_00", 1:m)
  cell_names <- paste0("cellline_", LETTERS[2:6], "A")
  values <- seq(0.1, 2.5, length.out = m * n)
  dt <- data.table::data.table(
    "Drug Name" = rep(drug_names, n),
    "Drug MOA" = rep(c(rep("moa_A", m - 1), "moa_B"), n),
    "Cell Line Name" = rep(cell_names, each = m),
    "Primary Tissue" = rep(c(rep("tissue_x", n - 1), "tissue_y"), each = m),
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
    "x" = values,
    "Concentration" = values,
    check.names = FALSE
  )
  list(drug_names = drug_names, cell_names = cell_names, dt = dt)
}