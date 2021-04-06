test_that("get_SE_experiment_metadata and set_SE_experiment_metadata work as expected", {
  exp_md <- list("Super" = "Star", "Serena" = "Williams")
  se <- SummarizedExperiment::SummarizedExperiment()

  oexp_md <- get_SE_experiment_metadata(se)
  expect_equal(oexp_md, NULL)
  se <- set_SE_experiment_metadata(se, exp_md)
  oexp_md <- get_SE_experiment_metadata(se)
  expect_equal(oexp_md, exp_md)
})


test_that("get_SE_keys and set_SE_keys work as expected", {
  keys <- list(Keys = list(Day0 = "TEST", Other = "STUFF"))
  se <- SummarizedExperiment::SummarizedExperiment(metadata = keys)
  nkeys <- get_SE_keys(se)
  expect_equal(nkeys$Day0, "TEST")
  expect_equal(nkeys$Other, "STUFF")

  # Test for all keys.
  keys2 <- list("test" = "NEW")
  se <- set_SE_keys(se, keys2)
  nkeys2 <- get_SE_keys(se)
  expect_equal(nkeys2$test, "NEW")
  expect_equal(length(nkeys2), 1)
})
