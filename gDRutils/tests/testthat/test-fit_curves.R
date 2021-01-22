#library(testthat); library(gDRutils)
source("setUp.R")
.round_params <- function(df) {
  df[] <- lapply(df, round, 4)
  df
}


test_that("fit_curves fails with expected errors", {
  expect_error(fit_curves(list()))

  # Log10 concentrations.
  df_resp2 <- df_resp
  df_resp2$Concentration <- df_resp2$Concentration * -1
  expect_error(fit_curves(df_resp2), regexp = "function accepts only unlogged concentrations")
})


test_that("warnings are thrown for duplicated concentrations", {
  df_resp3 <- df_resp
  df_resp3$Concentration[2:3] <- df_resp3$Concentration[1]
  expect_warning(fit_curves(df_resp3), regexp = "duplicate concentrations were found")
})


test_that("appropriate fit type is assigned for various use cases", {

  exp <- rbind(params, params_GR)
  rownames(exp) <- c("RV", "GR")

  # Test a 3P fit.
  df_result <- fit_curves(df_resp)
  expect_equal(.round_params(df_result[, names(params)]), exp, tolerance = 1e-5)

  obs_fit <- unique(df_result[, "fit_type"])
  expect_equal(obs_fit, "DRC3pHillFitModelFixS0")
  expect_equal(dim(df_result), c(2, 14))

  # Test a 4P fit (without the x_0 value).
  df_result <- fit_curves(df_resp, e_0 = NA, GR_0 = NA)
  expect_equal(.round_params(df_result[, names(params)]), exp, tolerance = 1e-5)
  obs_fit <- unique(df_result[, "fit_type"])
  expect_equal(obs_fit, "DRC4pHillFitModel")
  expect_equal(dim(df_result), c(2, 14))

  # Test for constant fit. 
  df_resp4 <- df_resp
  df_resp4$RelativeViability <- df_resp4$GRvalue <- 0.5

  expect_warning(df_result <- fit_curves(df_resp4), regexp = "overriding original x_0 argument") # Override. 
  obs_fit <- unique(df_result[ , "fit_type"])
  expect_equal(obs_fit, "DRCConstantFitResult")
  expect_equal(unname(unlist(df_result["RV", c("x_0", "x_inf", "x_mean", "x_AOC", "x_AOC_range")])), 
    rep(unique(df_resp4$RelativeViability), 5))
  expect_equal(dim(df_result), c(2, 14))

  # Test for all values above 0.5.
  # Test for all values below 0.5.
  # Test for a pushed constant fit. 
  # Test that force argument overrides as expected.

  # Test for too few points.
  df_result <- fit_curves(df_resp[3:5, ], n_point_cutoff = 4)
  obs_fit <- unique(df_result[, "fit_type"])
  expect_equal(obs_fit, "DRCTooFewPointsToFit")
  expect_equal(dim(df_result), c(2, 14))

  # TODO: Test for invalid fit. Maybe try a bunch of noise. 
#  expect_warning(df_result <- fit_curves(df_resp[3:5, ], n_point_cutoff = 1), regexp = "fitting failed with error")
#  
#  obs_fit <- unique(df_result[, "fit_type"])
#  expect_equal(obs_fit, "DRCInvalidFitResult")
#  expect_equal(dim(df_result), c(2, 14))
})
