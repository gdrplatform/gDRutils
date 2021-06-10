.round_params <- function(df) {
  df[] <- lapply(df, round, 4)
  df
}


test_that("fit_curves fails with expected errors", {
  expect_error(fit_curves(list()))

  # Log10 concentrations.
  df_resp2 <- df_resp
  df_resp2$Concentration <- df_resp2$Concentration * -1
  expect_error(fit_curves(df_resp2), reg = "function accepts only unlogged concentrations")

  # Invalid normalization_type.
  expect_error(fit_curves(df_resp, normalization_type = "BOGUS"), reg = "unknown curve type")
  expect_error(fit_curves(df_resp, normalization_type = c("GR", "BOGUS")), reg = "unknown curve type")
})


test_that("warnings are thrown for duplicated concentrations", {
  df_resp3 <- df_resp
  df_resp3$Concentration[2:3] <- df_resp3$Concentration[1]
  expect_warning(fit_curves(df_resp3), reg = "duplicate concentrations were found")
})


test_that("NA values are handled correctly", {
  df_resp_NA <- df_resp
  df_resp_NA[c(1, 4), c("RelativeViability", "GRvalue")] <- NA
  expect_error(fit_curves(df_resp_NA), NA)

  df_resp_NA2 <- df_resp_NA
  df_resp_NA2[6, "RelativeViability"] <- NA
  expect_error(fit_curves(df_resp_NA2), NA)

  df_resp_NA3 <- df_resp_NA
  df_resp_NA3[, "RelativeViability"] <- NA
  df_result_NA <- fit_curves(df_resp_NA3)
  expect_true(is.na(df_result_NA["RV", "xc50"]))
})


test_that("appropriate fit type is assigned for various use cases", {
  set.seed(1112020) # For reproducibility.

  # Test a 3P fit.
  ## Note that this should correspond to a cytotoxic response.
  df_result <- fit_curves(df_resp)
  expect_equal(.round_params(df_result[, names(params)]), expected, tolerance = 1e-5)

  obs_fit <- unique(df_result[, "fit_type"])
  expect_equal(obs_fit, "DRC3pHillFitModelFixS0")
  expect_equal(dim(df_result), expected_dims)

  # Test a 4P fit (without the x_0 value).
  df_result <- fit_curves(df_resp, e_0 = NA, GR_0 = NA)
  expect_equal(.round_params(df_result[, names(params)]), expected, tolerance = 1e-5)
  obs_fit <- unique(df_result[, "fit_type"])
  expect_equal(obs_fit, "DRC4pHillFitModel")
  expect_equal(dim(df_result), expected_dims)

  # Test for constant fit.
  df_resp4 <- df_resp
  df_resp4$RelativeViability <- df_resp4$GRvalue <- 0.5

  expect_warning(df_result <- fit_curves(df_resp4), reg = "overriding original x_0 argument") # Override.
  obs_fit <- unique(df_result[, "fit_type"])
  expect_equal(obs_fit, "DRCConstantFitResult")
  expect_equal(unname(unlist(df_result["RV", c("x_0", "x_inf", "x_mean", "x_AOC", "x_AOC_range")])),
    rep(unique(df_resp4$RelativeViability), 5))
  expect_equal(dim(df_result), expected_dims)

  ## Test for all values below 0.5.
  df_resp5 <- df_resp

  # Scale all readouts.
  max_rv <- max(df_resp$RelativeViability)
  min_rv <- min(df_resp$RelativeViability)

  max_gr <- max(df_resp$GRvalue)
  min_gr <- min(df_resp$GRvalue)

  df_resp5$RelativeViability <- (df_resp$RelativeViability / (max_rv - min_rv)) * ((max_rv - min_rv) / 2)
  df_resp5$GRvalue <- (df_resp$GRvalue / (max_gr - min_gr)) * ((max_gr - min_gr) / 2)
  df_result5 <- fit_curves(df_resp5)
  expect_equal(df_result5[, c("x_inf")], c(0, -1))
  obs_fit <- unique(df_result5[, "fit_type"])
  expect_equal(obs_fit, "DRC3pHillFitModelFixS0")

  # Test for all values above 0.5.
  ## Note that this corresponds to partial growth inhibition.
  df_result6 <- fit_curves(df_resp_above)
  expect_equal(df_result6[, c("x_inf")], c(0.5, 0.5), tolerance = 1e-5)
  expect_true(all(df_result6[, c("xc50")] > 500))
  obs_fit <- unique(df_result6[, "fit_type"])
  expect_equal(obs_fit, "DRC3pHillFitModelFixS0")

  # Test for a pushed constant fit by adding noise.
  df_resp7 <- df_resp_above
  noise <- sample(seq(-1, 1, 0.1), nrow(df_resp7))
  emax <- 0.8
  df_resp7$RelativeViability <- pmin(df_resp7$RelativeViability + noise, emax)
  df_resp7$GRvalue <- pmin(df_resp7$GRvalue + noise, emax)

  df_result7 <- fit_curves(df_resp7, force_fit = FALSE)
  expect_equal(unique(df_result7[, "fit_type"]), "DRCConstantFitResult")
  expect_equal(unique(unname(unlist(df_result7[, c("x_mean", "x_inf", "x_0")]))),
    0.70781, tolerance = 1e-5)

  # Test that force argument overrides as expected.
  df_result8 <- fit_curves(df_resp7, force_fit = TRUE)
  expect_equal(unique(df_result8[, "fit_type"]), "DRC3pHillFitModelFixS0")

  # Test that pcutoff argument works as expected.
  df_result9 <- fit_curves(df_resp7, force_fit = FALSE, pcutoff = 1)
  expect_equal(df_result9[, "fit_type"], c("DRC3pHillFitModelFixS0", "DRCConstantFitResult"))
  df_result10 <- fit_curves(df_resp7, force_fit = FALSE, pcutoff = 1.01) # Essentially equivalent to a 'force = TRUE'.
  expect_equal(unique(df_result10[, "fit_type"]), "DRC3pHillFitModelFixS0")

  # Test for GR values from 0-1.
  ## Note that this correspond to a fully cytostatic response (no growth).
  df_resp11 <- df_resp
  df_resp11$GRvalue <- df_resp11$RelativeViability

  df_result11 <- fit_curves(df_resp11)
  expect_equal(unique(df_result11[, "x_inf"]), c(0.1, 0.1), tolerance = 1e-5)

  # Test for too few points.
  df_result <- fit_curves(df_resp[3:5, ], n_point_cutoff = 4)
  obs_fit <- unique(df_result[, "fit_type"])
  expect_equal(obs_fit, "DRCTooFewPointsToFit")
  expect_equal(dim(df_result), expected_dims)

  #nolint start
    # TODO: Test for invalid fit. Maybe try a bunch of noise.
    #  expect_warning(df_result <- fit_curves(df_resp[3:5, ], n_point_cutoff = 1), reg = "fitting failed with error")
    #
    #  obs_fit <- unique(df_result[, "fit_type"])
    #  expect_equal(obs_fit, "DRCInvalidFitResult")
    #  expect_equal(dim(df_result), expected_dims)
  #nolint end
})


test_that("normalization_type can be specified", {
  GR_df_result <- fit_curves(df_resp, normalization_type = "GR")
  expect_equal(rownames(GR_df_result), "GR_gDR")
  expect_equal(.round_params(GR_df_result[, names(params)]), expected["GR_gDR", ], tolerance = 1e-5)

  RV_df_result <- fit_curves(df_resp, normalization_type = "RV")
  expect_equal(rownames(RV_df_result), "RV_gDR")
  expect_equal(.round_params(RV_df_result[, names(params)]), expected["RV_gDR", ], tolerance = 1e-5)
})


test_that(".estimate_xc50 works as expected", {
  expect_equal(gDRutils:::.estimate_xc50(c(NA, NA, NA, NA)), NA)
  expect_equal(gDRutils:::.estimate_xc50(c(0.6, 0.7, 0.8, 0.9)), Inf)
  expect_equal(gDRutils:::.estimate_xc50(c(0.1, 0.2, 0.3, 0.4)), -Inf)
  expect_equal(gDRutils:::.estimate_xc50(c(0.1, 0.2, 0.6, 0.7)), NA)
})


test_that("logistic_4parameters works as expected", {
    c <- 1
    Vinf <- 0.1
    V0 <- 1
    h <- 2
    EC50 <- 0.5

    # Non-numeric values cause an error.
    expect_error(logistic_4parameters(
        c = "non-numeric_entry",
        Vinf = Vinf,
        V0 = V0,
        EC50 = EC50,
        h = h
    ))

    # Normal fit.
    v <- logistic_4parameters(
        c = c,
        Vinf = Vinf,
        V0 = V0,
        EC50 = EC50,
        h = h
    )
    expect_equal(v, 0.28)

    # Flat fit.
    EC50 <- 0
    v <- logistic_4parameters(
        c = c,
        Vinf = Vinf,
        V0 = V0,
        EC50 = EC50,
        h = h
    )
    expect_equal(v, Vinf)
})


test_that(".setup_metric_output works as expected", {
  obs <- gDRutils:::.setup_metric_output()
  expect_true(is.list(obs))
  expect_equal(length(obs), 14)
})


test_that(".calculate_x_max works as expected", {
  n <- 4
  df <- data.frame(concs = c(0.01, 1, 0.03, 5), norm_values = seq(n)) 
  expect_equal(gDRutils:::.calculate_x_max(df), 2)

  df$norm_values <- rep(NA, n)
  expect_equal(gDRutils:::.calculate_x_max(df), NA)
})
