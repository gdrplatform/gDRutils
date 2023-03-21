listMAE <- lapply(list.files(system.file(package = "gDRtestData", "testdata"),
                             "final", full.names = TRUE)[1:2], readRDS)
listSE <- lapply(listMAE, function(x) x[[2]])
names(listSE) <- c("combo1", "combo2")


test_that("merge_assay works as expected", {
  normalizedMerged <- merge_assay(listSE, "Normalized")
  checkmate::expect_list(normalizedMerged)
  testthat::expect_true(all(c("DT", "BM") == names(normalizedMerged)))
  checkmate::expect_class(normalizedMerged[[1]], "data.table")
  checkmate::expect_class(normalizedMerged[[2]], "BumpyDataFrameMatrix")
})

test_that("merge_metadata and identify_unique_se_metadata_fields work as expected", {
  metadata_fields <- identify_unique_se_metadata_fields(listSE)
  mergedMetadata <- merge_metadata(listSE, metadata_fields)
  expect_identical(names(mergedMetadata), metadata_fields)
  expect_identical(names(mergedMetadata$experiment_metadata), names(listSE))
  
  listSE2 <- listSE
  newMetaName <- "dummy_meta"
  S4Vectors::metadata(listSE2$combo1)[[newMetaName]] <- list()
  metadata_fields2 <- identify_unique_se_metadata_fields(listSE2)
  expect_true(newMetaName %in% metadata_fields2)
  mergedMetadata2 <- merge_metadata(listSE2, metadata_fields2)
  expect_true(newMetaName %in% names(mergedMetadata2))
})

test_that("merge_SE works as expected", {
  mergedSE <- purrr::quietly(merge_SE)(listSE, discard_keys = c("normalization_type", "fit_source"))
  expect_length(mergedSE$warnings, 2)
  checkmate::expect_class(mergedSE$result, "SummarizedExperiment")
  S4Vectors::metadata(mergedSE$result)[["df_raw_data"]] <- list(NULL)
  validate_SE(mergedSE$result)
  additional_col_name <- "QCS"
  mergedSE2 <- purrr::quietly(merge_SE)(listSE, additional_col_name, discard_keys = c("normalization_type", "fit_source"))
  expect_length(mergedSE2$warnings, 2)
  assayNormalized <- convert_se_assay_to_dt(mergedSE2$result, "Metrics") 
  expect_true(additional_col_name %in% names(assayNormalized))
  expect_identical(unique(assayNormalized[[additional_col_name]]), names(listSE))
  })
