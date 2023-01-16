context("standardize_MAE")


mapping_vector <- as.character(seq_len(5))
names(mapping_vector) <- names(iris)

test_that("rename_bumpy works as expected",  {
  bumpy <- BumpyMatrix::splitAsBumpyMatrix(iris, iris$Species, iris$Sepal.Length)
  bumpy_renamed <- rename_bumpy(bumpy, mapping_vector)
  expect_equal(BumpyMatrix::commonColnames(bumpy_renamed), unname(mapping_vector))
})

test_that("rename_DFrame works as expected",  {
  dframe <- S4Vectors::DataFrame(iris)
  dframe_renamed <- rename_DFrame(dframe, mapping_vector)
  expect_equal(names(dframe_renamed), unname(mapping_vector))
})


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

