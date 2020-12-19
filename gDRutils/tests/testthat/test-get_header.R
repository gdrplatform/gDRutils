#library(testthat); library(gDRutils)

test_that("get_header works", {
  print("test get_header is running")
  expect_error(get_header("BOGUS"))
  expect_equal(get_header("manifest"), c("Barcode", "Template", "Duration"))
  expect_equal(length(get_header()), 13)
})
