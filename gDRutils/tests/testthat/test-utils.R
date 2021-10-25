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
