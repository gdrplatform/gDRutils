library(testthat)

context("assay_to_dt")

test_that("assay_to_dt works as expected", {
  SE <- readRDS(system.file(package ="gDRutils", "inst", "testdata", "exemplarySE.rds"))
  normalized <- assay_to_dt(SE, "Normalized")
  data.table::setcolorder(normalized, sort(names(normalized)))
  data.table::setorder(normalized)
  normalizedRef <- data.table::as.data.table(readRDS(system.file(package ="gDRutils", "inst", "testdata", "normalizedExemplary.rds")))
  data.table::setcolorder(normalizedRef, sort(names(normalizedRef)))
  data.table::setorder(normalizedRef)
  testthat::expect_equal(normalized, normalizedRef)
})

