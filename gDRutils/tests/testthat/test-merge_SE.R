test_that("merge_assay works as expected", {
  listSE <- lapply(list.files(system.file(package = "gDRtestData", "testdata"), "final", full.names = TRUE)[1:2], readRDS)
  names(listSE) <- c("combo1", "combo2")
  normalizedMerged <- merge_assay(listSE, "Normalized")
  checkmate::expect_list(normalizedMerged, names = c("DT", "BM"))
  checkmate::expect_class(normalizedMerged[[1]], "data.table")
  checkmate::expect_class(normalizedMerged[[2]], "BumpyDataFrameMatrix")
})

test_that("merge_SE works as expected", {
  listSE <- lapply(list.files(system.file(package = "gDRtestData", "testdata"), "final", full.names = TRUE)[1:2], readRDS)
  names(listSE) <- c("combo1", "combo2")
  mergedSE <- merge_SE(listSE)
  S4Vectors::metadata(mergedSE)
  checkmate::expect_class(mergedSE, "SummarizedExperiment")
  validate_SE(mergedSE)
  additional_col_name <- "QCS"
  mergedSE2 <- merge_SE(listSE, additional_col_name)
  assayNormalized <- convert_se_assay_to_dt(mergedSE2, "Metrics") 
  expect_true(additional_col_name %in% names(assayNormalized))
  expect_identical(unique(assayNormalized[[additional_col_name]]), names(listSE))
  })