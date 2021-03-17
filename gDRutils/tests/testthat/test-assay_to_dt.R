library(testthat)
library(gDRutils)

test_that("convert_se_assay_to_dt works as expected", {
  m <- 20
  n <- 10
  rnames <- LETTERS[1:m]
  cnames <- letters[1:n]

  # Normal matrix.
  ref_gr_value <-  matrix(runif(m * n), nrow = m, ncol = n, dimnames = list(rnames, cnames))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(RefGRvalue = ref_gr_value),
                                                   rowData = S4Vectors::DataFrame(rnames),
                                                   colData = S4Vectors::DataFrame(cnames))

  dt <- convert_se_assay_to_dt(se = se, assay_name = "RefGRvalue", include_metadata = FALSE)
  expect_equal(dt$RefGRvalue, as.vector(ref_gr_value))
  expect_equal(dim(dt), c(200, 3))

  dt <- convert_se_assay_to_dt(se = se, assay_name = "RefGRvalue", include_metadata = TRUE)
  expect_equal(dt[order(dt$cId), "RefGRvalue"][[1]], as.vector(ref_gr_value))
  expect_equal(dim(dt), c(200, 5))
  expect_equal(dt$rnames, as.character(dt$rId))
  expect_equal(dt$cnames, as.character(dt$cId))

  # BumpyDataFrameMatrix.
  df <- S4Vectors::DataFrame(r = rep(rnames, n), c = rep(cnames, m), values = runif(m * n))
  norm <- BumpyMatrix::splitAsBumpyMatrix(df, row = df$r, column = df$c)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(norm = norm),
                                                   rowData = S4Vectors::DataFrame(rnames),
                                                   colData = S4Vectors::DataFrame(cnames))
  dt <- convert_se_assay_to_dt(se = se, assay_name = "norm", include_metadata = FALSE)
  merged <- merge(df, S4Vectors::DataFrame(dt[, c("rId", "cId", "values")]))
  expect_equal(merged$r, merged$rId)
  expect_equal(merged$c, merged$cId)

  dt <- convert_se_assay_to_dt(se = se, assay_name = "norm", include_metadata = TRUE)
  merged <- merge(df, S4Vectors::DataFrame(dt[, c("rnames", "cnames", "values")]))
  expect_equal(merged$r, merged$rnames)
  expect_equal(merged$c, merged$cnames)
})


test_that("assay_to_dt works as expected", {
  # DataFrameMatrix.
  SE <- readRDS(system.file(package = "gDRutils", "testdata", "exemplarySE.rds"))
  normalized <- assay_to_dt(SE, "Normalized")
  data.table::setcolorder(normalized, sort(names(normalized)))
  data.table::setorder(normalized)
  normalizedRef <- 
    data.table::as.data.table(
      readRDS(system.file(package = "gDRutils", "testdata", "normalizedExemplary.rds"))
    )
  data.table::setcolorder(normalizedRef, sort(names(normalizedRef)))
  data.table::setorder(normalizedRef)
  testthat::expect_equal(normalized, normalizedRef)
})
