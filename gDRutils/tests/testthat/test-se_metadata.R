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


test_that("get_SE_fit_parameters and set_SE_fit_parameters work as expected", {
  params <- list(n_point_cutoff = 10,
                 range_conc = c(1, 100),
                 force_fit = TRUE,
                 pcutoff = 1,
                 cap = 0.2)
  se <- SummarizedExperiment::SummarizedExperiment(metadata = list())
  fit_params <- get_SE_fit_parameters(se)
  
  expect_equal(fit_params, NULL)

  se <- set_SE_fit_parameters(se, params)
  expect_equal(get_SE_fit_parameters(se), params)
})


test_that("get_SE_identifiers works as expected", {
  exp <- list("drug" = "drug", "celllinename" = "CellLineName")
  se <- SummarizedExperiment::SummarizedExperiment(metadata = list(identifiers = exp))

  obs <- get_SE_identifiers(se)
  expect_equal(obs, exp)

  obs <- get_SE_identifiers(se, "celllinename")
  expect_equal(obs, exp[["celllinename"]])

  # Identifier does not exist on the SummarizedExperiment. 
  # Used mainly for backwards compatibility purposes.
  obs <- get_SE_identifiers(se, "masked_tag")
  expect_equal(obs, "masked")

  # Invalid identifier.
  expect_error(get_SE_identifiers(se, "INVALID"))
})
