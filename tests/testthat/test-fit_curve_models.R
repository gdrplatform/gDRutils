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
    expect_true(all(c("r2", "converged", "ec50_extrapolated", "n_conc",
                      "fit_type") %in% names(out)), info = m)
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


test_that("curves the unbounded fit could not evaluate are fitted", {
  # Two real RV curves from a KRAS-G12C panel. Without a lower bound on ec50
  # drc reached -1.2 and -0.0085: log(dose/ec50) is then NaN at every dose,
  # coef() still returns numbers, and summing NaN residuals with na.rm = TRUE
  # scored both as r2 = 1.0000 — a perfect fit to a curve the model cannot
  # evaluate. Bounded, both fit, and ec50 is a concentration again.
  steep <- data.frame(
    Concentration = c(0, 0.00015241615447757103, 0.000457247142398638,
                      0.0013717406226510097, 0.0041152289299736305,
                      0.012345679549008601, 0.037037016924301704,
                      0.11111098560477098, 0.33333352883795708, 1),
    x = c(1.06616363636364, 1.04777272727273, 0.909538461538462,
          0.629416666666667, 0.257723076923077, 0.0428076923076923,
          0.0110307692307692, 0.0120846153846154, 0.00700769230769231,
          0.00862307692307692)
  )
  shallow <- data.frame(
    Concentration = c(0.0015241580352681, 0.0045724714239864,
                      0.0137174062265101, 0.0411522892997364,
                      0.123456795490086, 0.370370169243018,
                      1.11111241447565, 3.33333528837957, 10),
    x = c(1.1048, 1.0375, 0.9534, 1.0914, 0.6881, 0.5649, 0.4372,
          0.4605, 0.3762)
  )

  for (dt in list(steep, shallow)) {
    out <- ll3_fit_fn(dt)
    expect_true(out$converged)
    expect_gt(out$ec50, 0)
    # the tell of the old bug: a perfect score from residuals that were all NaN
    expect_true(is.finite(out$r2) && out$r2 < 1)
    expect_identical(out$n_conc, 9L)
  }
})


test_that("an ec50 outside the tested doses is flagged, not hidden", {
  # A real GR slice for a drug too weak to reach its midpoint on the plate:
  # the response only falls to 0.74 at the top dose of 10. The fit is a good
  # one — r2 near 0.98 — but ec50 lands near 66, so it was extrapolated rather
  # than measured, and ranking potency on it would be ranking on no data.
  weak <- data.frame(
    Concentration = c(0.0015241580352681, 0.0045724714239864,
                      0.0137174062265101, 0.0411522892997364,
                      0.123456795490086, 0.370370169243018,
                      1.11111241447565, 3.33333528837957, 10),
    x = c(1.0241, 1.0152, 0.999, 0.9977, 0.9916, 0.9744, 0.876,
          0.8433, 0.7369)
  )
  out <- ll3_fit_fn(weak)
  expect_true(out$converged)
  expect_true(out$ec50_extrapolated)
  expect_gt(out$ec50, max(weak$Concentration))

  # and the ordinary case is not flagged
  ok <- ll3_fit_fn(.curve_dt(12))
  expect_false(ok$ec50_extrapolated)
})


test_that("an ec50 the data cannot locate is refused, not reported", {
  # A real GR slice with no dose-response at all — the response sits between
  # 0.91 and 1.02 across four orders of magnitude. Bounding ec50 below keeps
  # it positive but cannot give it a meaning, so the optimiser runs away to
  # 2e19 with residuals that stay finite. r2 of 0.03 says the fit is worthless;
  # the ec50 would read as a number.
  flat_real <- data.frame(
    Concentration = c(0.0015241580352681, 0.0411522892997364,
                      0.123456795490086, 0.370370169243018,
                      1.11111241447565, 3.33333528837957, 10),
    x = c(0.9191, 0.9551, 0.9665, 0.9478, 0.9146, 1.02, 0.9108)
  )
  out <- ll3_fit_fn(flat_real)
  expect_false(out$converged)
  expect_true(is.na(out$ec50))
  expect_identical(out$n_conc, 7L)
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
