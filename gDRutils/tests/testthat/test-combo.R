library(testthat)
context("combo-related functions")
test_se <-
  readRDS(gDRtestData::get_test_dataset_paths()["finalSE_combo_with_metrics"])
test_l <- convert_combo_data_to_dt(test_se)

test_that("convert_combo_data_to_dt",  {
  res_l <- convert_combo_data_to_dt(test_se)

  # expected assays converted
  exp_as <- as.character(get_combo_assay_names())
  expect_identical(sort(names(res_l)), sort(exp_as))

  # check content of data.table
  expect_true(nrow(res_l[[1]]) > 1 && ncol(res_l[[1]]) > 1)
  exp_idfs <- get_prettified_identifiers(c("drugname", "drugname2", "cellline"), simplify = FALSE)
  expect_true(all(exp_idfs %in% colnames(res_l[[1]])))

  # errors
  expect_error(
    convert_combo_data_to_dt(data.frame(a = 1)),
    "Assertion on 'se' failed: Must inherit from class 'SummarizedExperiment', but has class 'data.frame'.",
    fixed = TRUE
  )
  err_msg <-
    "Assertion on 'dummy, Smooth_Matrix' failed. Must be element(s) of {'SmoothMatrix, BlissExcess, "
  err_msg2 <-
    "HSAExcess, HSAScore, BlissScore, CIScore_50, CIScore_80, isobolograms'} set."
  expect_error(
    convert_combo_data_to_dt(test_se, c_assays = c("dummy", "Smooth_Matrix", "HSAExcess")),
    sprintf("%s%s", err_msg, err_msg2),
    fixed = TRUE
  )
  expect_error(
    convert_combo_data_to_dt(test_se, normalization_type = 1),
    "Assertion on 'normalization_type' failed: Must be of type 'character', not 'double'.",
    fixed = TRUE
  )
  expect_error(
    convert_combo_data_to_dt(test_se, prettify = "true"),
    "Assertion on 'prettify' failed: Must be of type 'logical flag', not 'character'.",
    fixed = TRUE
  )
})
