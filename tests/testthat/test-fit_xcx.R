# A clean decreasing sigmoid with a low plateau, so most levels are reachable.
.xcx_dt <- function(x_inf = 0.05) {
  conc <- 10^seq(-3, 2, length.out = 10)
  data.table::data.table(
    Concentration = conc,
    x = x_inf + (1 - x_inf) / (1 + (conc / 0.5)^1.5),
    normalization_type = "RV"
  )
}


test_that("every requested level becomes a column named for it", {
  out <- fit_drug_response_metrics_xcx(.xcx_dt(), levels = c(75, 90))
  expect_true(all(c("xc75", "xc90") %in% names(out)))
  # the fit's own columns survive untouched
  expect_true(all(c("xc50", "ec50", "h", "x_inf", "x_0", "fit_type") %in% names(out)))
  expect_true(all(lengths(out) == 1L))
})


test_that("the number is inhibition, not the response level", {
  # xc90 is where the response reaches 0.1, so it must sit ABOVE xc50 on a
  # decreasing curve. Under the opposite reading it would sit below, which is
  # the whole reason this is asserted rather than assumed.
  out <- fit_drug_response_metrics_xcx(.xcx_dt(), levels = c(75, 90))
  expect_gt(out$xc75, out$xc50)
  expect_gt(out$xc90, out$xc75)
})


test_that("level 50 returns the fit's own xc50, not a second column", {
  # Both invert the same curve, so a second xc50 could only disagree or
  # duplicate. Asking for it is satisfied by the value already there.
  out <- fit_drug_response_metrics_xcx(.xcx_dt(), levels = c(50, 90))
  expect_identical(sum(names(out) == "xc50"), 1L)
  expect_true(is.finite(out$xc50))
})


test_that("the helper reproduces the fit's own xc50", {
  # The decisive check: .calculate_xc50() and this helper both invert the same
  # model, so they must agree. Asserted on a fit whose xc50 is finite and
  # inside the tested range, where no capping is in play.
  fit <- fit_drug_response_metrics(.xcx_dt())
  skip_if_not(is.finite(fit$xc50), "no finite xc50 to compare against")
  bare <- fit[setdiff(names(fit), "xc50")]
  again <- gDRutils:::.xc_at_levels(bare, levels = 50)
  expect_equal(again$xc50, fit$xc50, tolerance = 1e-8)
})


test_that("a level the curve never reaches is Inf, not a number", {
  # Plateau at 0.4, so 90% inhibition is unreachable and 50% is not.
  out <- fit_drug_response_metrics_xcx(.xcx_dt(x_inf = 0.4), levels = c(90))
  expect_false(is.finite(out$xc90))
  expect_true(is.na(out$xc90) || out$xc90 == Inf)
})


test_that("a fit that failed gives NA levels rather than raising", {
  too_few <- data.table::data.table(
    Concentration = c(0.1, 1),
    x = c(0.9, 0.5),
    normalization_type = "RV"
  )
  expect_error(fit_drug_response_metrics_xcx(too_few, levels = c(90)), NA)
  out <- fit_drug_response_metrics_xcx(too_few, levels = c(90))
  expect_true("xc90" %in% names(out))
  expect_false(is.finite(out$xc90))
})


test_that("levels are validated, because they become column names", {
  dt <- .xcx_dt()
  expect_error(fit_drug_response_metrics_xcx(dt, levels = 0))
  expect_error(fit_drug_response_metrics_xcx(dt, levels = 100))
  expect_error(fit_drug_response_metrics_xcx(dt, levels = 37.5))
  expect_error(fit_drug_response_metrics_xcx(dt, levels = c(90, 90)))
  expect_error(fit_drug_response_metrics_xcx(dt, levels = numeric(0)))
  expect_error(fit_drug_response_metrics_xcx(dt, levels = NA_integer_))
})


test_that("it carries a description for the provenance", {
  info <- attr(fit_drug_response_metrics_xcx, "fit_fn_info")
  expect_false(is.null(info))
  expect_true(nzchar(info$model))
  expect_true(nzchar(info$description))
})


test_that("it is registered for single-agent data", {
  expect_true("fit_drug_response_metrics_xcx" %in%
                get_fit_fns(data_type = "single-agent")$name)
  expect_identical(get_fit_fn("fit_drug_response_metrics_xcx"),
                   fit_drug_response_metrics_xcx)
})
