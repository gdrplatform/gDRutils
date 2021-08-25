test_that("get_header works", {
  reset_env_identifiers()

  expect_error(get_header("BOGUS"))
  expect_equal(get_header("manifest"), c("Barcode", "Template", "Duration"))
  expect_equal(length(get_header()), 12)

  set_env_identifier("duration", "TEST_DURATION")
  expect_equal(get_header("manifest"), c("Barcode", "Template", "TEST_DURATION"))

  reset_env_identifiers()
  expect_equal(get_header("manifest"), c("Barcode", "Template", "Duration"))
})
