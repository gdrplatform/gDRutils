test_that(".clean_key_inputs works as expected", {
  keys <- LETTERS[1:5]
  cols <- LETTERS[1:3]
  expect_warning(out <- gDRutils:::.clean_key_inputs(keys, cols))
  expect_equal(out, cols)
})


test_that("assert_equal_input_len works as expected", {
  ec50 <- 0.5
  x_0 <- 1
  x_inf <- 0.1
  h <- 2
  efficacy <- 0.6
  expect_equal(gDRutils:::assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h), NULL)

  # Bad lengths.
  ec50 <- c(0.5, 0.5)
  expect_error(gDRutils:::assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h))

  # Length 1 fit parameters.
  ec50 <- 0.5
  efficacy <- c(0.6, 0.7, 0.8)
  expect_equal(gDRutils:::assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h), NULL)

  # Length 1 outlier.
  ec50 <- c(0.5, 0.6)
  x_0 <- c(1, 0.9)
  x_inf <- c(0.1, 0.15)
  h <- c(2, 2)
  efficacy <- 0.6
  expect_equal(gDRutils:::assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h), NULL)
})

test_that("assert_choices",  {
  ### expected values
  expect_null(assert_choices(letters[1], letters))
  expect_null(assert_choices(letters[1:2], letters))
  expect_null(assert_choices(1:5, 1:10))
  expect_null(assert_choices(1, 1:10))
  
  ### errors
  err_msg <- sprintf("Assertion on '%s' failed.", letters[8])
  expect_error(assert_choices(letters[5:8], letters[1:7]), err_msg)
  err_msg <-
    sprintf("Assertion on '%s, %s' failed.", letters[7], letters[8])
  expect_error(assert_choices(letters[5:8], letters[1:6]), err_msg)
})

test_that("MAEpply works as expectd", {
  maeReal <- readRDS(system.file("testdata", "finalMAE_combo_2dose_nonoise2.RDS", package = "gDRtestData"))
  list1 <- MAEpply(maeReal, SummarizedExperiment::assayNames)
  expect_length(list1, 2)
  expect_true(inherits(list1, "list"))
  list2 <- MAEpply(maeReal, SummarizedExperiment::assayNames, unify = TRUE)
  expect_length(list1, 5)
  expect_true(inherits(list2, "character"))
})


