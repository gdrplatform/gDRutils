library(testthat)

context("assay_to_df")

test_that("assay_to_df works as expected", {
  SE <- readRDS(system.file(package ="gDRutils", "inst", "testdata", "exemplarySE.rds"))
  normalized <- assay_to_df(SE, "Normalized")
  normalizedRef <- readRDS(system.file(package ="gDRutils", "inst", "testdata", "normalizedExemplary.rds"))
  testthat::expect_identical(normalized, normalizedRef)
})
