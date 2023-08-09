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

test_that("standardize_se works as expected",  {
  se_original <- get_synthetic_data("finalMAE_combo_matrix_small")[[1]]
  se <- se_original
  se@metadata$identifiers$drug <- "druuug"
  se@metadata$identifiers$concentration2 <- "dose 2"
  rowData(se) <- rename_DFrame(rowData(se), c("Gnumber" = "druuug"))
  assay(se, "RawTreated") <- rename_bumpy(assay(se, "RawTreated"), c("Concentration_2" = "dose 2"))
  se_standardized <- standardize_se(se)
  expect_equal(convert_se_assay_to_dt(se_standardized, "RawTreated"),
               convert_se_assay_to_dt(se_original, "RawTreated"))
})


test_that("standardize_MAE works as expected",  {
  mae_original <- get_synthetic_data("finalMAE_combo_matrix_small")
  mae <- mae_original
  mae[[1]]@metadata$identifiers$drug <- "druuug"
  mae[[2]]@metadata$identifiers$drug <- "druuug"
  mae[[1]]@metadata$identifiers$concentration2 <- "dose 2"
  rowData(mae[[1]]) <- rename_DFrame(rowData(mae[[1]]), c("Gnumber" = "druuug"))
  rowData(mae[[2]]) <- rename_DFrame(rowData(mae[[2]]), c("Gnumber" = "druuug"))
  assay(mae[[1]], "RawTreated") <- rename_bumpy(assay(mae[[1]], "RawTreated"), c("Concentration_2" = "dose 2"))
  mae_standardized <- standardize_mae(mae)
  expect_equal(convert_mae_assay_to_dt(mae_standardized, "RawTreated"),
               convert_mae_assay_to_dt(mae_original, "RawTreated"))
})

test_that("standardize_MAE works with polymapped identifiers",  {
  mae_original <- get_synthetic_data("finalMAE_combo_matrix_small")
  mae <- mae_original
  mae[[1]]@metadata$identifiers$drug <- c("druuug", "Drug")
  mae[[2]]@metadata$identifiers$drug <- c("druuug", "Drug")
  mae[[1]]@metadata$identifiers$concentration2 <- "dose 2"
  rowData(mae[[1]]) <- rename_DFrame(rowData(mae[[1]]), c("Gnumber" = "druuug"))
  rowData(mae[[2]]) <- rename_DFrame(rowData(mae[[2]]), c("Gnumber" = "druuug"))
  assay(mae[[1]], "RawTreated") <- rename_bumpy(assay(mae[[1]], "RawTreated"), c("Concentration_2" = "dose 2"))
  mae_standardized <- standardize_mae(mae)
  expect_equal(convert_mae_assay_to_dt(mae_standardized, "RawTreated"),
               convert_mae_assay_to_dt(mae_original, "RawTreated"))
})

test_that("colData/rowData refinement functions work as expected",  {
  mae <- get_synthetic_data("finalMAE_combo_matrix_small")
  expect_true(inherits(refine_coldata(SummarizedExperiment::colData(mae[[1]]), mae[[1]]),
               "DataFrame"))
  expect_true(inherits(refine_rowdata(SummarizedExperiment::rowData(mae[[1]]), mae[[1]]),
               "DataFrame"))
  
  expect_error(refine_coldata(mae, mae), "Assertion on 'se' failed:")
  expect_error(refine_rowdata(mae, mae), "Assertion on 'se' failed:")
  
})
