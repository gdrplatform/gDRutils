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


test_that("standardize_MAE works as expected",  {
  mae_original <- readRDS(system.file(package = "gDRtestData", "testdata", "finalMAE_combo_matrix_small.RDS"))
  mae <- mae_original
  mae[[1]]@metadata$identifiers$drug <- "druuug"
  mae[[2]]@metadata$identifiers$drug <- "druuug"
  mae[[1]]@metadata$identifiers$concentration2 <- "dose 2"
  rowData(mae[[1]]) <- rename_DFrame(rowData(mae[[1]]), c("Gnumber" = "druuug"))
  rowData(mae[[2]]) <- rename_DFrame(rowData(mae[[2]]), c("Gnumber" = "druuug"))
  assay(mae[[1]], "RawTreated") <- rename_bumpy(assay(mae[[1]], "RawTreated"), c("Concentration_2" = "dose 2"))
  mae_standardized <- standardize_MAE(mae)
  expect_equal(get_MAE_identifiers(mae_standardized),
               get_MAE_identifiers(mae_original))
  expect_equal(convert_mae_assay_to_dt(mae_standardized, "RawTreated"),
               convert_mae_assay_to_dt(mae_original, "RawTreated"))
})
