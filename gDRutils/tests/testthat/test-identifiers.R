
test_that("get_identifier and set_identifier work", {
  reset_identifiers()

  expect_error(get_identifier("BOGUS"))

  expect_equal(get_identifier("duration"), "Duration")
  
  expect_equal(get_identifier(c("duration", "barcode")), c("Duration", "Barcode"))
  
  set_identifier("cellline", "my_personal_cell_line_identifiers")
  expect_equal(get_identifier("cellline"), "my_personal_cell_line_identifiers")
  vals <- gDRutils::IDENTIFIERS_LIST
  vals[["cellline"]] <- "my_personal_cell_line_identifiers"

  expect_equal(get_identifier(), vals)
})


test_that("reset_identifiers works", {
  reset_identifiers()

  d <- get_identifier()
  expect_true(length(d) != 0L)
  set_identifier("duration", "TEST_DURATION")
  expect_equal(get_identifier("duration"), "TEST_DURATION")

  ilist <- reset_identifiers()
  expect_equal(ilist, NULL)
  expect_equal(get_identifier("duration"), "Duration")
})
