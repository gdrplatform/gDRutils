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
  expect_warning( # Bartek Czech confirmed this warning is expected
    se <- set_SE_keys(se, keys2),
    "overwriting existing metadata entry: 'Keys'"
  )
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
  exp <- list("drug" = "drug", "celllinename" = "CellLineName", "buggy_idfs" = "test", "masked_tag" = "masked")
  se <- SummarizedExperiment::SummarizedExperiment(metadata = list(identifiers = exp))
  expect_equal(get_SE_identifiers(se), exp, simplify = FALSE)
  expect_error(get_SE_identifiers(se, "buggy_idfs", simplify = TRUE), 
               "Assertion on 'id_type' failed: Must be element of set")
  expect_error(get_SE_identifiers(se, "INVALID", simplify = TRUE), 
               "Assertion on 'id_type' failed: Must be element of set")

  # Identifier does not exist on the SummarizedExperiment,
  # so get it from the environment.
  obs <- get_SE_identifiers(se, "masked_tag", simplify = TRUE)
  expect_equal(obs, "masked")

  # Set identifiers.
  expect_warning(
    se <- set_SE_identifiers(se, list()),
    "overwriting existing metadata entry: 'identifiers'"
  )
  obs2 <- get_SE_identifiers(se, simplify = TRUE)
  expect_equal(obs2, list())

  # Multiple identifiers.
  exp <- list("drug_name" = "Drugs", "cellline_name" = "Cells", "duration" = "Duration")
  se <- SummarizedExperiment::SummarizedExperiment(metadata = list(identifiers = exp))
  expect_equal(get_SE_identifiers(se, c("drug_name", "duration"), simplify = FALSE), 
               list(drug_name = "Drugs", duration = "Duration")) # Env and se identifiers.
  expect_equal(get_SE_identifiers(se, c("cellline_name", "drug_name"), simplify = FALSE), 
               list(cellline_name = "Cells", drug_name = "Drugs")) # Order.
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

test_that("get_SE_experiment_raw_data and set_SE_experiment_raw_data work as expected", {
  se <- SummarizedExperiment::SummarizedExperiment(metadata = list())
  experiment_raw_data <- get_SE_experiment_raw_data(se)
  expect_equal(experiment_raw_data, NULL)
  raw_data <- data.frame(a = 1:3, b = 3:5)
  se <- set_SE_experiment_raw_data(se, raw_data)
  expect_equal(get_SE_experiment_raw_data(se), raw_data)
})

