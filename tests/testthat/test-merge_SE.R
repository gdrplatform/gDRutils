listMAE <- lapply(list.files(system.file(package = "gDRtestData", "testdata"),
                             "final", full.names = TRUE)[1:2], qs::qread)
listSE <- lapply(listMAE, function(x) x[[2]])
names(listSE) <- c("combo1", "combo2")

listMAE2 <- lapply(list.files(system.file(package = "gDRtestData", "testdata"),
                             "final", full.names = TRUE)[1:2], qs::qread)
listSE2 <- lapply(listMAE, function(x) x[[1]])
names(listSE2) <- c("combo1", "combo2")

listSE3 <- c(listSE[1], listSE2[2])

listSE4 <- c(listSE[1], ligand = qs::qread(
  list.files(system.file(package = "gDRtestData", "testdata"),
             "Ligand", full.names = TRUE))[[1]])


test_that("merge_assay works as expected", {
  normalizedMerged <- merge_assay(listSE, "Normalized")
  checkmate::expect_list(normalizedMerged)
  testthat::expect_true(all(c("DT", "BM") == names(normalizedMerged)))
  checkmate::expect_data_table(normalizedMerged[[1]])
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
  set_env_identifier("cellline", "CELL_LINE_ID")
  mergedSE <- purrr::quietly(merge_SE)(listSE)
  checkmate::expect_class(mergedSE$result, "SummarizedExperiment")
  S4Vectors::metadata(mergedSE$result)[["df_raw_data"]] <- list(NULL)
  validate_SE(mergedSE$result)
  additional_col_name <- "QCS"
  mergedSE2 <- purrr::quietly(merge_SE)(listSE, additional_col_name)
  assayNormalized <- convert_se_assay_to_dt(mergedSE2$result, "Metrics") 
  expect_true(additional_col_name %in% names(assayNormalized))
  expect_identical(unique(assayNormalized[[additional_col_name]]), names(listSE))
  expect_identical(SummarizedExperiment::assayNames(listSE[[1]]),
                   SummarizedExperiment::assayNames(mergedSE[[1]]))
  reset_env_identifiers()
  })


test_that("merge_SE works as expected with combo matrix data", {
  mergedSE <- purrr::quietly(merge_SE)(listSE2)
  checkmate::expect_class(mergedSE$result, "SummarizedExperiment")
  S4Vectors::metadata(mergedSE$result)[["df_raw_data"]] <- list(NULL)
  validate_SE(mergedSE$result)
  additional_col_name <- "QCS"
  mergedSE2 <- purrr::quietly(merge_SE)(listSE2, additional_col_name)
  assayNormalized <- convert_se_assay_to_dt(mergedSE2$result, "Metrics") 
  expect_true(additional_col_name %in% names(assayNormalized))
  expect_identical(unique(assayNormalized[[additional_col_name]]), names(listSE))
  expect_identical(SummarizedExperiment::assayNames(listSE2[[1]]),
                   SummarizedExperiment::assayNames(mergedSE[[1]]))
})

test_that("merge_SE works as expected with mixed data types", {
  mergedSE <- purrr::quietly(merge_SE)(listSE3)
  checkmate::expect_class(mergedSE$result, "SummarizedExperiment")
  S4Vectors::metadata(mergedSE$result)[["df_raw_data"]] <- list(NULL)
  validate_SE(mergedSE$result)
})


test_that("merge_SE works with data with additional perturbations", {
  mergedSE <- purrr::quietly(merge_SE)(listSE4)
  checkmate::expect_class(mergedSE$result, "SummarizedExperiment")
  validate_SE(mergedSE$result)
  expect_equal(dim(mergedSE$result), c(10, 5))
})
