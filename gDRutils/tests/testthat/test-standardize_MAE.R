context("standardize_MAE")

test_that("colData/rowData refinement functions work as expected",  {
  mae <-
    readRDS(
      system.file(package = "gDRtestData", "testdata", "finalMAE_combo_matrix_small.RDS")
    )
  expect_true(inherits(refine_coldata(SummarizedExperiment::colData(mae[[1]]), mae[[1]]),
               "DataFrame"))
  expect_true(inherits(refine_rowdata(SummarizedExperiment::rowData(mae[[1]]), mae[[1]]),
               "DataFrame"))
  
  expect_error(refine_coldata(mae, mae), "Assertion on 'se' failed:")
  expect_error(refine_rowdata(mae, mae), "Assertion on 'se' failed:")
  
})
