test_that("get_env_identifiers and set_env_identifier work", {
  reset_env_identifiers()

  expect_error(get_env_identifiers("BOGUS", simplify = TRUE))
  expect_equal(get_env_identifiers("duration", simplify = TRUE), "Duration")
  expect_equal(get_env_identifiers(c("duration", "barcode"), simplify = FALSE), 
               list(duration = "Duration", barcode = "Barcode"))
  
  set_env_identifier("cellline", "my_personal_cell_line_identifiers")
  expect_equal(get_env_identifiers("cellline", simplify = TRUE), "my_personal_cell_line_identifiers")
  vals <- gDRutils:::IDENTIFIERS_LIST
  vals[["cellline"]] <- "my_personal_cell_line_identifiers"

  expect_equal(get_env_identifiers(simplify = TRUE), vals)
})


test_that("reset_env_identifiers works", {
  reset_env_identifiers()

  d <- get_env_identifiers(simplify = TRUE)
  set_env_identifier("duration", "TEST_DURATION")
  expect_equal(get_env_identifiers("duration", simplify = TRUE), "TEST_DURATION")

  ilist <- reset_env_identifiers()
  expect_equal(ilist, NULL)
  expect_equal(get_env_identifiers("duration", simplify = TRUE), "Duration")
})


test_that("get_prettified_identifiers works as expected", {
  obs <- get_prettified_identifiers(c("drug_name", "cellline_name"), simplify = FALSE)
  expect_equal(obs, c("Drug Name", "Cell Line Name"))

  obs <- get_prettified_identifiers(simplify = TRUE)
  expect_true(is(obs, "list"))
  expect_true(length(obs) > 1L)
})

test_that("get_env_identifier works with untreated_tag", {
  obs <- get_env_identifiers("untreated_tag", simplify = TRUE)
  expect_length(obs, 2)
  expect_equal(obs, c("untreated", "vehicle"))
})

test_that("get_SE_identifier works with untreated_tag", {
  se <- SummarizedExperiment()
  obs <- get_SE_identifiers(se, "untreated_tag", simplify = TRUE)
  expect_length(obs, 2)
  expect_equal(obs, c("untreated", "vehicle"))
})

test_that("get_expect_one_identifiers works as expected", {
  expect_true(all(get_expect_one_identifiers() %in% names(get_env_identifiers())))
  expect_true(length(get_expect_one_identifiers()) > 1L)
})

test_that("get_required_identifiers works as expected", {
  expect_true(all(get_required_identifiers() %in% names(get_env_identifiers())))
  expect_true(length(get_required_identifiers()) > 1L)
  expect_true(all(get_required_identifiers() %in% get_expect_one_identifiers()))
})
