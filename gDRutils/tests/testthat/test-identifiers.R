test_that("get_env_identifiers and set_env_identifier work", {
  reset_env_identifiers()

  expect_error(get_env_identifiers("BOGUS"))
  expect_equal(get_env_identifiers("duration"), "Duration")
  expect_equal(get_env_identifiers(c("duration", "barcode")), c("Duration", "Barcode"))
  
  set_env_identifier("cellline", "my_personal_cell_line_identifiers")
  expect_equal(get_env_identifiers("cellline"), "my_personal_cell_line_identifiers")
  vals <- gDRutils:::IDENTIFIERS_LIST
  vals[["cellline"]] <- "my_personal_cell_line_identifiers"

  expect_equal(get_env_identifiers(), vals)
})


test_that("reset_env_identifiers works", {
  reset_env_identifiers()

  d <- get_env_identifiers()
  expect_true(length(d) != 0L)
  set_env_identifier("duration", "TEST_DURATION")
  expect_equal(get_env_identifiers("duration"), "TEST_DURATION")

  ilist <- reset_env_identifiers()
  expect_equal(ilist, NULL)
  expect_equal(get_env_identifiers("duration"), "Duration")
})


test_that("support deprecated get_identifiers", {
  obs <- get_identifier("drugname")
  expect_equal(obs, "DrugName")
})
