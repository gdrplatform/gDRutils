#library(testthat); library(gDRutils)

context("test functions relating to getting and resetting headers")

test_that("get_header works", {
  expect_error(get_header("BOGUS"))
  expect_equal(get_header("manifest"), c("Barcode", "Template", "Duration"))
  expect_equal(length(get_header()), 13)
})


test_that("reset_headers works", {
  set_identifier("duration", "TEST_DURATION")
  expect_equal(get_header("manifest"), c("Barcode", "Template", "TEST_DURATION"))

  hlist <- reset_headers()
  expect_equal(get_header("manifest"), c("Barcode", "Template", "Duration"))
  expect_equal(hlist, NULL)
})
