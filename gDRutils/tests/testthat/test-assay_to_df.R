library(testthat)

context("assay_to_df")

test_that("assay_to_df works as expected", {
  SE <- readRDS(system.file(package ="gDRutils", "inst", "testdata", "exemplarySE.rds"))
  normalized <- assay_to_df(SE, "Normalized")
  normalized <- normalized[, order(colnames(normalized))]
  normalizedRef <- data.table::as.data.table(readRDS(system.file(package ="gDRutils", "inst", "testdata", "normalizedExemplary.rds")))
  normalizedRef <- normalizedRef[, order(colnames(normalizedRef))]
  testthat::expect_identical(sort(normalized), sort(normalizedRef))
})
