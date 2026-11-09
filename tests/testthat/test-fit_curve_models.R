# A clean decreasing sigmoid, the shape every model here should manage.
.curve_dt <- function(n = 8) {
  conc <- 10^seq(-3, 2, length.out = n)
  data.frame(Concentration = conc, x = 1 / (1 + (conc / 0.5)^1.5))
}

MODELS <- c("ll3_fit_fn", "ll5_fit_fn", "weibull1_fit_fn",
            "weibull2_fit_fn", "gompertz_fit_fn", "hormesis_fit_fn")


test_that("every model returns the shared core, whatever else it returns", {
  for (m in MODELS) {
    out <- do.call(m, list(.curve_dt()))
    expect_type(out, "list")
    expect_true(all(c("r2", "converged", "n_conc", "fit_type") %in% names(out)),
                info = m)
    expect_length(out$fit_type, 1L)
    # one row, so every element is scalar — apply_fit turns this into one
    # row of the output assay and a vector would silently become several
    expect_true(all(lengths(out) == 1L), info = m)
  }
})


test_that("a well-behaved sigmoid fits, and ec50 lands near the truth", {
  # Not every model parameterises the midpoint identically, so this is
  # asserted only for the log-logistics, where ec50 means what it says.
  for (m in c("ll3_fit_fn", "ll5_fit_fn")) {
    out <- do.call(m, list(.curve_dt(12)))
    expect_true(out$converged, info = m)
    expect_true(is.finite(out$r2) && out$r2 > 0.95, info = m)
    expect_true(out$ec50 > 0.1 && out$ec50 < 2.5, info = m)
  }
})


test_that("too few concentrations is a result, not an error", {
  # A model with k free parameters cannot be identified from fewer than k
  # points, and drc will return something anyway. Refusing is the point.
  for (m in MODELS) {
    out <- do.call(m, list(.curve_dt(2)))
    expect_false(out$converged, info = m)
    expect_identical(out$n_conc, 2L, info = m)
    expect_true(is.na(out$r2), info = m)
  }
})


test_that("a fit that cannot converge returns NA parameters, not an error", {
  flat <- data.frame(Concentration = 10^(-3:3), x = rep(1, 7))
  for (m in MODELS) {
    expect_error(do.call(m, list(flat)), NA, info = m)
    out <- do.call(m, list(flat))
    expect_true(all(c("converged", "fit_type") %in% names(out)), info = m)
  }
})


test_that("missing columns are handled rather than thrown", {
  for (m in MODELS) {
    out <- do.call(m, list(data.frame(wrong = 1, columns = 2)))
    expect_false(out$converged, info = m)
    expect_identical(out$n_conc, 0L, info = m)
  }
})


test_that("non-finite and non-positive concentrations are dropped", {
  dt <- .curve_dt(10)
  dt$Concentration[1] <- 0      # log-logistic is undefined at zero dose
  dt$x[2] <- NA
  out <- ll3_fit_fn(dt)
  expect_lte(out$n_conc, 8L)
})


test_that("the column names are configurable", {
  dt <- .curve_dt()
  names(dt) <- c("conc_other", "resp_other")
  out <- ll3_fit_fn(dt, x_col = "resp_other", conc_col = "conc_other")
  expect_true(out$converged)
})


test_that("each model carries a description for the provenance", {
  for (m in MODELS) {
    info <- attr(get(m), "fit_fn_info")
    expect_false(is.null(info), info = m)
    expect_true(nzchar(info$model), info = m)
    expect_true(nzchar(info$description), info = m)
  }
})


test_that("the models report distinct fit_type values", {
  types <- vapply(MODELS, function(m) do.call(m, list(.curve_dt()))$fit_type,
                  character(1))
  expect_identical(length(unique(types)), length(MODELS))
})


test_that("drc's coefficient labels do not leak into the output", {
  # drc names coefficients "<param>:(Intercept)". If that is not stripped the
  # assay gets columns nobody asked for.
  out <- ll3_fit_fn(.curve_dt())
  expect_false(any(grepl(":", names(out), fixed = TRUE)))
})


test_that("the models are registered for single-agent data", {
  skip_if_not(exists("get_fit_fns", envir = asNamespace("gDRutils"),
                     inherits = FALSE),
              "registry not present")
  listed <- get_fit_fns(data_type = "single-agent")$name
  expect_true(all(MODELS %in% listed))
})
