dataSEpath <- lapply(list.files(system.file(package = "gDRtestData", "testdata"),
                                "finalSE", full.names = TRUE)[1:2], readRDS)

names(dataSEpath) <- c("combo1", "combo2")


test_that("merge_assay works as expected", {
  additional_col_name <- "synthetic"
  mergedNorm <- gDRutils::merge_assay(dataSEpath, assay_name = "Metrics", additional_col_name = additional_col_name)
  expect_length(mergedNorm, 2)
  expect_named(mergedNorm, c("DT", "BM"))
  expect_true(additional_col_name %in% names(mergedNorm$DT))
})
  