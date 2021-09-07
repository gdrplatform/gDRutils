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


test_that("get_prettified_identifiers works as expected", {
  obs <- get_prettified_identifiers(c("drugname", "cellline_name"))
  expect_equal(obs, c("Drug", "Cell Line"))

  obs <- get_prettified_identifiers()
  expect_true(is(obs, "list"))
  expect_true(length(obs) > 1L)
})

test_that("get_identifier works with untreated_tag", {
  obs <- get_env_identifiers("untreated_tag")
  expect_length(obs, 2)
  expect_equal(obs, c("untreated", "vehicle"))
})
