library(testthat)
context("combo-related functions")

test_that("convert_combo_data_to_dt", {
  test_mae <- get_synthetic_data("finalMAE_combo_matrix_small")
  res_l <- convert_combo_data_to_dt(test_mae[[1]])

  # expected assays converted
  exp_as <- as.character(get_combo_assay_names())
  expect_identical(sort(names(res_l)), sort(exp_as))

  # check content of data.table
  expect_true(nrow(res_l[[1]]) > 1 && ncol(res_l[[1]]) > 1)
  exp_idfs <- get_prettified_identifiers(c("drug_name", "drug_name2", "cellline"), simplify = FALSE)
  expect_true(all(exp_idfs %in% colnames(res_l[[1]])))

  # errors
  expect_error(
    convert_combo_data_to_dt(data.table::data.table(a = 1)),
    paste("Assertion on 'se' failed: Must inherit from class 'SummarizedExperiment',",
          "but has classes 'data.table','data.frame'."), 
    fixed = TRUE
  )
  err_msg <-
    "Assertion on 'dummy, Smooth_Matrix' failed. Must be element(s) of {'SmoothMatrix, BlissExcess, "
  err_msg2 <-
    "HSAExcess, HSAScore, BlissScore, CIScore_50, CIScore_80, isobolograms'} set."
  expect_error(
    convert_combo_data_to_dt(test_mae[[1]], c_assays = c("dummy", "Smooth_Matrix", "HSAExcess")),
    sprintf("%s%s", err_msg, err_msg2),
    fixed = TRUE
  )
  expect_error(
    convert_combo_data_to_dt(test_mae[[1]], normalization_type = 1),
    "Assertion on 'normalization_type' failed: Must be of type 'character', not 'double'.",
    fixed = TRUE
  )
  expect_error(
    convert_combo_data_to_dt(test_mae[[1]], prettify = "true"),
    "Assertion on 'prettify' failed: Must be of type 'logical flag', not 'character'.",
    fixed = TRUE
  )
})

test_that("get_iso_colors",  {
  ### expected values
  gic <- get_iso_colors()
  expect_true(length(gic) > 2)
  expect_identical("character", class(gic))
  gic2 <- get_iso_colors(formals(get_iso_colors)[[1]][[3]])
  expect_true(any(gic != gic2))
  expect_identical(length(gic), length(gic2))
  
  ### errors
  expect_error(get_iso_colors("inv_param"), "'arg' should be one of ")
})

test_that("get_combo_col_settings",  {
  ### expected values
  gcan <- names(get_combo_assay_names()[1])
  gcc <-
    get_combo_col_settings(g_metric = "GR", assay_type = gcan)
  expect_true(inherits(gcc, "list"))
  expect_identical(sort(names(gcc)), c("breaks", "colors", "limits"))
  
  ### errors
  err_msg <- "Assertion on 'assay_type' failed: "
  expect_error(get_combo_col_settings("GR", 8), err_msg)
  err_msg <- "Assertion on 'g_metric' failed: "
  expect_error(get_combo_col_settings("grvalue", 8), err_msg)
})

test_that("shorten_normalization_type_name", {
  ### expected values
  expect_identical("GR", shorten_normalization_type_name("GRvalue"))
  
  ### errors
  err_msg <- "Assertion on 'x' failed: Must be element of set"
  expect_error(shorten_normalization_type_name("invalid"), err_msg)
})
