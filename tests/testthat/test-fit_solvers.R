context("Hill fit solvers")

.solver_dt <- function(x = c(0.98, 0.95, 0.88, 0.70, 0.45, 0.25, 0.14, 0.10), norm_type = "RV") {
  data.table::data.table(
    Concentration = 10 ^ seq(-3, 1, length.out = length(x)),
    x = x,
    x_std = 0.02,
    normalization_type = norm_type
  )
}

# The same bounded Hill fit by stats::nls(algorithm = "port") — a second optimiser,
# used to check that the solver really is a parameter.
.solver_nls <- function(conc, x, model, control = list()) {
  fixed_x_0 <- model$fixed$x_0
  f <- if (is.null(fixed_x_0)) {
    x ~ x_inf + (x_0 - x_inf) / (1 + (conc / ec50) ^ h)
  } else {
    x ~ x_inf + (fixed_x_0 - x_inf) / (1 + (conc / ec50) ^ h)
  }
  fit <- tryCatch(
    stats::nls(f, data = data.frame(conc = conc, x = x), start = as.list(model$start),
               lower = model$lower, upper = model$upper, algorithm = "port",
               control = utils::modifyList(list(warnOnly = TRUE), control)),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(params = c(h = NA, x_inf = NA, x_0 = NA, ec50 = NA), rss = NA,
                converged = FALSE, message = conditionMessage(fit)))
  }
  p <- stats::coef(fit)
  list(params = c(h = p[["h"]], x_inf = p[["x_inf"]],
                  x_0 = if (is.null(fixed_x_0)) p[["x_0"]] else fixed_x_0, ec50 = p[["ec50"]]),
       rss = sum(stats::residuals(fit) ^ 2), converged = isTRUE(fit$convInfo$isConv),
       message = NA_character_)
}


test_that("solver_drc is the default and is registered", {
  expect_identical(formals(fit_drug_response_metrics)$solver, quote(solver_drc))
  expect_true("drc" %in% names(get_solvers()))
  expect_identical(get_solvers()[["drc"]], solver_drc)
})


test_that("solver_drc returns the parameters fit_drug_response_metrics reports", {
  dt <- .solver_dt()
  fit <- fit_drug_response_metrics(dt)
  model <- gDRutils:::.hill_model(
    1, start = c(2, 0.4, stats::median(dt$Concentration)),
    lower = c(0.1, 0, min(dt$Concentration) / 10),
    upper = c(5, 1, max(dt$Concentration) * 10)
  )
  sol <- solver_drc(dt$Concentration, dt$x, model)
  expect_true(sol$converged)
  expect_equal(unname(sol$params[c("h", "x_inf", "ec50")]), c(fit$h, fit$x_inf, fit$ec50), tolerance = 1e-8)
  expect_equal(sol$params[["x_0"]], 1)
  expect_equal(sol$rss, fit$rss, tolerance = 1e-8)
})


test_that("the predictions come from the Hill closed form, as drc's would", {
  dt <- .solver_dt()
  fit <- fit_drug_response_metrics(dt)
  grid <- 10 ^ seq(log10(min(dt$Concentration)), log10(max(dt$Concentration)), length.out = 101L)
  drc_fit <- drc::drm(x ~ conc, data = data.frame(conc = dt$Concentration, x = dt$x),
                      fct = drc::LL.3u(upper = 1, names = c("h", "x_inf", "ec50")),
                      start = c(2, 0.4, stats::median(dt$Concentration)),
                      lowerl = c(0.1, 0, min(dt$Concentration) / 10),
                      upperl = c(5, 1, max(dt$Concentration) * 10))
  expect_equal(fit$x_mean, mean(stats::predict(drc_fit, data.frame(conc = grid))), tolerance = 1e-8)
})


test_that("another solver can replace drc and finds the same curve", {
  for (norm_type in c("RV", "GR")) {
    dt <- .solver_dt(norm_type = norm_type)
    ref <- fit_drug_response_metrics(dt)
    got <- fit_drug_response_metrics(dt, solver = .solver_nls)
    expect_equal(got$fit_type, ref$fit_type, info = norm_type)
    expect_equal(c(got$h, got$x_inf, got$ec50, got$xc50), c(ref$h, ref$x_inf, ref$ec50, ref$xc50),
                 tolerance = 1e-3, info = norm_type)
  }
  # the 4-parameter model too
  ref4 <- fit_drug_response_metrics_4p(.solver_dt())
  got4 <- fit_drug_response_metrics_4p(.solver_dt(), solver = .solver_nls)
  expect_equal(got4$fit_type, "DRC4pHillFitModel")
  expect_equal(c(got4$x_0, got4$ec50), c(ref4$x_0, ref4$ec50), tolerance = 1e-3)
})


test_that("a solver that fails or throws gives an invalid fit, not an error", {
  dt <- .solver_dt()
  not_converged <- function(conc, x, model, control = list()) {
    list(params = c(h = NA, x_inf = NA, x_0 = NA, ec50 = NA), rss = NA,
         converged = FALSE, message = "no luck")
  }
  throws <- function(conc, x, model, control = list()) stop("boom")
  expect_equal(fit_drug_response_metrics(dt, solver = not_converged)$fit_type, "DRCInvalidFitResult")
  expect_warning(fit <- fit_drug_response_metrics(dt, solver = throws), "boom")
  expect_equal(fit$fit_type, "DRCInvalidFitResult")
  # a solver that only reports converged = FALSE does not warn
  expect_silent(fit_drug_response_metrics(dt, solver = not_converged))
})


test_that("solver_control reaches the solver", {
  seen <- NULL
  spy <- function(conc, x, model, control = list()) {
    seen <<- control
    solver_drc(conc, x, model, control)
  }
  fit_drug_response_metrics(.solver_dt(), solver = spy, solver_control = list(relTol = 1e-8))
  expect_equal(seen, list(relTol = 1e-8))
})


test_that("x_0 fixes the upper asymptote, as the combination co-treatment fits need", {
  dt <- .solver_dt(x = c(0.78, 0.75, 0.70, 0.55, 0.36, 0.20, 0.12, 0.08))
  fit <- fit_drug_response_metrics(dt, x_0 = 0.8)
  expect_equal(fit$x_0, 0.8)
  expect_equal(fit$fit_type, "DRC3pHillFitModelFixS0")
})


test_that("a solver passed through apply_fit is recorded and replayed", {
  mae <- get_synthetic_data("finalMAE_small.qs2")
  se <- mae[["single-agent"]][1:2, 1:2]
  se <- apply_fit(se, fit_fn = fit_drug_response_metrics, data_type = "single-agent",
                  output_assay = "nls_fit", fit_source = "nls",
                  fit_args = list(solver = .solver_nls))
  rec <- get_fit_provenance(se, "nls_fit", "nls")
  expect_true("solver" %in% names(rec$fn_args))
  chk <- check_reproducibility(se, "nls_fit", "nls", allow_source = TRUE)
  expect_true(all(chk$reproducible), info = paste(chk$detail, collapse = "; "))
})


test_that("register_solver adds to the menu", {
  register_solver("nls_port_test", .solver_nls, "nls port, for the tests")
  on.exit(rm("nls_port_test", envir = gDRutils:::.solver_registry), add = TRUE)
  expect_identical(get_solvers()[["nls_port_test"]], .solver_nls)
  register_solver("never_available", .solver_nls, available = function() FALSE)
  on.exit(rm("never_available", envir = gDRutils:::.solver_registry), add = TRUE)
  expect_false("never_available" %in% names(get_solvers()))
  expect_true("never_available" %in% names(get_solvers(available_only = FALSE)))

  info <- get_solver_info()
  expect_identical(names(info), c("name", "description", "available"))
  expect_identical(sort(info$name), sort(names(get_solvers(available_only = FALSE))))
  expect_false(info$available[info$name == "never_available"])
  expect_identical(info$description[info$name == "nls_port_test"], "nls port, for the tests")
  expect_false("never_available" %in% get_solver_info(available_only = TRUE)$name)
})


test_that("solver_drc reports a silent drc failure as not converged", {
  dt <- .solver_dt()
  nan_start <- function(conc, x, model, control = list()) {
    model$start[] <- NaN
    solver_drc(conc, x, model, control)
  }
  model <- gDRutils:::.hill_model(1, start = c(2, 0.4, 0.1), lower = c(0.1, 0, 1e-4), upper = c(5, 1, 100))
  sol <- suppressWarnings(nan_start(dt$Concentration, dt$x, model))
  expect_false(sol$converged)
  expect_true(all(is.na(sol$params)))
  fit <- suppressWarnings(fit_drug_response_metrics(dt, solver = nan_start))
  expect_equal(fit$fit_type, "DRCInvalidFitResult")
})


test_that("a solver must take control, and its errors are not silent", {
  dt <- .solver_dt()
  no_control <- function(conc, x, model) solver_drc(conc, x, model)
  expect_error(fit_drug_response_metrics(dt, solver = no_control), "control")
  expect_error(register_solver("no_control", no_control), "control")
  broken <- function(conc, x, model, control = list()) not_a_solver_function(conc)
  expect_warning(fit <- fit_drug_response_metrics(dt, solver = broken), "not_a_solver_function")
  expect_equal(fit$fit_type, "DRCInvalidFitResult")
})
