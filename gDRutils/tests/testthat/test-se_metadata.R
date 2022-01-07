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


test_that("get_SE_identifiers and set_SE_identifiers works as expected", {
  exp <- list("drug" = "gDrug", "cellline_name" = "gCell")
  se <- SummarizedExperiment::SummarizedExperiment(metadata = list(identifiers = exp))

  # No identifier passed.
  obs <- get_SE_identifiers(se, simplify = TRUE)
  expect_equal(obs, exp)
  
  # Single identifier.
  obs <- get_SE_identifiers(se, "cellline_name", simplify = TRUE)
  expect_equal(obs, exp[["cellline_name"]])
  
  # Invalid identifier.
  exp <- list("drug" = "drug", "celllinename" = "CellLineName", "buggy_idfs" = "test")
  se <- SummarizedExperiment::SummarizedExperiment(metadata = list(identifiers = exp))
  expect_equal(get_SE_identifiers(se), exp, simplify = FALSE)
  expect_error(get_SE_identifiers(se, "buggy_idfs", simplify = TRUE), 
               "Assertion on 'id_type' failed: Must be element of set")
  expect_error(get_SE_identifiers(se, "INVALID", simplify = TRUE), 
               "Assertion on 'id_type' failed: Must be element of set")

  # Identifier does not exist on the SummarizedExperiment,
  # so get it from the environment.
  expect_warning(obs <- get_SE_identifiers(se, "masked_tag", simplify = TRUE), 
    regexp = "'se' was passed, but identifier 'masked_tag' not found on se's identifiers")
  expect_equal(obs, "masked")

  # Set identifiers.
  se <- set_SE_identifiers(se, list())
  obs2 <- get_SE_identifiers(se, simplify = TRUE)
  expect_equal(obs2, list())

  # Multiple identifiers.
  exp <- list("drugname" = "Drugs", "cellline_name" = "Cells")
  se <- SummarizedExperiment::SummarizedExperiment(metadata = list(identifiers = exp))
  expect_equal(get_SE_identifiers(se, c("drugname", "duration"), simplify = FALSE), 
               list(drugname = "Drugs", duration = "Duration")) # Env and se identifiers.
  expect_equal(get_SE_identifiers(se, c("cellline_name", "drugname"), simplify = FALSE), 
               list(cellline_name = "Cells", drugname = "Drugs")) # Order.
})

test_that("get_SE_processing_metadata and set_SE_processing_metadata work as expected", {
  params <- list(date_processed = Sys.Date(),
                 session_info = sessionInfo())
  se <- SummarizedExperiment::SummarizedExperiment(metadata = list())
  processing_metadata <- get_SE_processing_metadata(se)
  
  expect_equal(processing_metadata, NULL)
  
  se <- set_SE_processing_metadata(se, params)
  expect_equal(get_SE_processing_metadata(se), params)
})

