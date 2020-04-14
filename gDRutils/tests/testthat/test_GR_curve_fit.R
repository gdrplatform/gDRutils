library(testthat)

context("ICGRfits")

test_that("ICGRfits fails with expected errors", {
  expect_error(ICGRfits(list()),
  regexp = "inherits\\(df_, \"data.frame\"\\) is not TRUE")
})
