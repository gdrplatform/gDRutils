test_that("prettify_flat_metrics works as expected", {
  x <- c("CellLineName", "Tissue",
         "Primary Tissue",
         "GR_gDR_x_mean", "GR_gDR_xc50", 
         "RV_GDS_x_mean", 
         "Concentration_2", "Gnumber_2", "Drug_3",
         "E_0", "GR_gDR_x_AOC_range"
  )
  
  y <- c("Gnumber", "Gnumber_2", "MyDrug", "MyDrug_2")
  
  obs <- prettify_flat_metrics(x, human_readable = FALSE)
  exp <- c("CellLineName", "Tissue",
           "Primary Tissue",
           "GR_mean", "GR50", 
           "GDS_RV_mean", 
           "Concentration_2", "Gnumber_2", "Drug_3",
           "E_0", "GR_AOC_range")
  expect_equal(obs, exp)
  
  # Human readable names work.
  obs <- prettify_flat_metrics(x, human_readable = TRUE)
  exp <- c("Cell Line Name", "Tissue",
           "Primary Tissue",
           "GR Mean", "GR50", 
           "RV Mean (GDS)",
           "Concentration 2", "Gnumber 2", "Drug 3",
           "E0", "GR AOC within set range")
  expect_equal(obs, exp)
  obs <- prettify_flat_metrics(y, human_readable = TRUE)
  exp <- c("Gnumber", "Gnumber 2", "My Drug", "My Drug 2")
  expect_equal(obs, exp)
})

test_that(".convert_norm_specific_metrics works as expected", {
  
  idfs <- get_env_identifiers(k = NULL, simplify = TRUE)
  norm_type <- c("GR", "RV")
  expect_equal(.convert_norm_specific_metrics(idfs, norm_type), idfs)
  expect_equal(.convert_norm_specific_metrics(idfs, c("AB", "BC")), idfs)
  
  col_name_1 <- c(
    "rId", "cId", "GR_gDR_x_mean", "GR_gDR_x_AOC", "GR_gDR_x_AOC_range", "GR_gDR_xc50",
    "GR_gDR_x_max", "GR_gDR_ec50", "GR_gDR_x_inf", "GR_gDR_x_0", "GR_gDR_h", "GR_gDR_r2",
    "GR_gDR_x_sd_avg", "GR_gDR_fit_type", "maxlog10Concentration", "N_conc", "cotrt_value", "source",
    "Gnumber", "DrugName", "drug_moa", "Gnumber_2", "DrugName_2", "drug_moa_2",
    "Duration", "clid", "CellLineName", "Tissue", "parental_identifier", "subtype",
    "ReferenceDivisionTime", "RV_gDR_x_mean", "RV_gDR_x_AOC", "RV_gDR_x_AOC_range", "RV_gDR_xc50", "RV_gDR_x_max",
    "RV_gDR_ec50", "RV_gDR_x_inf", "RV_gDR_x_0", "RV_gDR_h", "RV_gDR_r2", "RV_gDR_x_sd_avg",
    "RV_gDR_fit_type"
  )
  col_name_1_exp <- c(
    "rId", "cId", "_gDR_GR_mean", "_gDR_GR_AOC", "_gDR_GR_AOC_range", "_gDR_GR50",
    "_gDR_GR_max", "_gDR_GEC50", "_gDR_GR_inf", "_gDR_GR_0", "_gDR_h_GR", "_gDR_GR_r2",
    "_gDR_GR_sd_avg", "_gDR_fit_type_GR", "maxlog10Concentration", "N_conc", "cotrt_value", "source",
    "Gnumber", "DrugName", "drug_moa", "Gnumber_2", "DrugName_2", "drug_moa_2",
    "Duration", "clid", "CellLineName", "Tissue", "parental_identifier", "subtype",
    "ReferenceDivisionTime", "_gDR_RV_mean", "_gDR_RV_AOC", "_gDR_RV_AOC_range", "_gDR_IC50", "_gDR_E_max",
    "_gDR_EC50", "_gDR_E_inf", "_gDR_E_0", "_gDR_h_RV", "_gDR_RV_r2", "_gDR_RV_sd_avg",
    "_gDR_fit_type_RV"
  )
  expect_equal(.convert_norm_specific_metrics(col_name_1, norm_type), col_name_1_exp)
  
  
  col_name_2 <- c(
    "rId", "cId", "Concentration", "Gnumber", "DrugName", "drug_moa",
    "Duration", "clid", "CellLineName", "Tissue", "parental_identifier", "subtype",
    "ReferenceDivisionTime", "RelativeViability", "GRvalue", "std_RelativeViability", "std_GRvalue"
  )
  expect_equal(.convert_norm_specific_metrics(col_name_2, norm_type), col_name_2)
  expect_equal(.convert_norm_specific_metrics(col_name_2, c("AB", "BC")), col_name_2)
  
  col_name_combo_1 <- c(
    "rId", "cId", "Gnumber", "DrugName", "drug_moa", "Gnumber_2",
    "DrugName_2", "drug_moa_2", "Duration", "clid", "CellLineName", "Tissue",
    "parental_identifier", "subtype", "ReferenceDivisionTime", "hsa_score_GR", "hsa_score_RV", "bliss_score_GR",
    "bliss_score_RV", "CIScore_50_GR", "CIScore_50_RV", "CIScore_80_GR", "CIScore_80_RV"
  )
  expect_equal(.convert_norm_specific_metrics(col_name_combo_1, norm_type), col_name_combo_1)
  expect_equal(.convert_norm_specific_metrics(col_name_combo_1, c("AB", "BC")), col_name_combo_1)
  
  col_name_combo_2 <- c(
    "rId", "cId", "Concentration", "Concentration_2", "Gnumber", "DrugName",
    "drug_moa", "Gnumber_2", "DrugName_2", "drug_moa_2", "Duration", "clid",
    "CellLineName", "Tissue", "parental_identifier", "subtype", "ReferenceDivisionTime", "smooth_GR",
    "smooth_RV", "hsa_excess_GR", "hsa_excess_RV", "bliss_excess_GR", "bliss_excess_RV"
  )
  expect_equal(.convert_norm_specific_metrics(col_name_combo_2, norm_type), col_name_combo_2)
  expect_equal(.convert_norm_specific_metrics(col_name_combo_2, c("AB", "BC")), col_name_combo_2)
  
  col_name_3 <- c(
    "GR_gDR_x_max", "GR_gDR_ec50", "smooth_GR", "bliss_excess_GR",
    "RV_gDR_x_max", "RV_gDR_ec50", "smooth_RV", "bliss_excess_RV"
  )
  col_name_3_exp <- c(
    "GR_gDR_x_max", "GR_gDR_ec50", "smooth_GR", "bliss_excess_GR",
    "_gDR_E_max", "_gDR_EC50", "smooth_RV", "bliss_excess_RV"
  )
  expect_equal(.convert_norm_specific_metrics(col_name_3, c("AB", "RV")), col_name_3_exp)
  
})

test_that(".prettify_metadata_columns works as expected", {
  col_name <- c(
    "_gDR_EC50", "_gDR_E_inf", "_gDR_E_0", "_gDR_h_RV",
    "ReferenceDivisionTime", "RelativeViability", "GRvalue", "std_RelativeViability", "std_GRvalue",
    "smooth_GR", "hsa_excess_GR", "bliss_excess_GR"
  )
  col_name_exp <- c(
    "EC50", "E Inf", "E 0", "h RV",
    "Reference Division Time", "Relative Viability", "GRvalue", "Std Relative Viability", "Std GRvalue", 
    "Smooth GR", "HSA Excess GR", "Bliss Excess GR"
  )
  expect_equal(.prettify_metadata_columns(col_name), col_name_exp)
  
})

