#library(testthat); library(gDRutils)

test_that("get_header works", {
  reset_identifiers()
  reset_headers()

  expect_error(get_header("BOGUS"))
  expect_equal(get_header("manifest"), c("Barcode", "Template", "Duration"))
  expect_equal(length(get_header()), 13)
})


test_that("reset_headers works", {
  reset_identifiers()
  reset_headers()

  set_identifier("duration", "TEST_DURATION")
  expect_equal(get_header("manifest"), c("Barcode", "Template", "TEST_DURATION"))

  reset_identifiers()
  hlist <- reset_headers()
  expect_equal(get_header("manifest"), c("Barcode", "Template", "Duration"))
  expect_equal(hlist, NULL)
})
