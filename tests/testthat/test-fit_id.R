suppressWarnings(
  se_id <- qs2::qs_read(
    system.file("testdata/finalMAE_small.qs2", package = "gDRtestData")
  )[["single-agent"]]
)
S4Vectors::metadata(se_id)$provenance <- NULL

# a function standing in for user code: see test-fit_provenance.R
.as_user <- function(fn, values = list()) {
  environment(fn) <- list2env(values, parent = globalenv())
  fn
}

.fit_hill <- function(se, fit_source, ..., slicing_values = NULL, fit_args = list()) {
  suppressWarnings(apply_fit(
    se, fit_drug_response_metrics, "single-agent",
    output_assay = "hill", fit_source = fit_source,
    slicing_values = slicing_values, fit_args = fit_args, ...
  ))
}

.hill_ids <- function(se, norm) {
  dt <- gDRutils:::.assay_dt(se, "hill")
  unique(as.character(dt[["fit_id"]][as.character(dt[["normalization_type"]]) == norm]))
}


test_that("every row written by apply_fit carries the fit_id of its recipe", {
  se <- .fit_hill(se_id, "loose", fit_args = list(pcutoff = 0.1))
  ids <- unique(as.character(gDRutils:::.assay_dt(se, "hill")[["fit_id"]]))

  expect_length(ids, 1L)
  expect_match(ids, "^[0-9a-f]{12}$")
  expect_identical(names(get_SE_provenance(se)$fits$hill), ids)

  recipe <- get_fit_recipe(se, "hill", ids)
  expect_identical(recipe$fit_id, ids)
  expect_identical(recipe$fit_source, "loose")
  expect_identical(get_fit_provenance(se, "hill", "loose")$fit_id, ids)
  expect_identical(get_SE_provenance(se)$schema_version, 2L)
})


test_that("the fit_id names the label and what a fit does, not its slices or data", {
  loose <- .fit_hill(se_id, "loose", fit_args = list(pcutoff = 0.1))
  relabelled <- .fit_hill(se_id, "renamed", fit_args = list(pcutoff = 0.1))
  one_slice <- .fit_hill(se_id, "loose", slicing_values = "GR", fit_args = list(pcutoff = 0.1))
  stricter <- .fit_hill(se_id, "loose", fit_args = list(pcutoff = 0.5))

  id_of <- function(se) names(get_SE_provenance(se)$fits$hill)
  # the same label, code and settings are the same fit, whatever slices were fitted
  expect_identical(id_of(one_slice), id_of(loose))
  # another label or other settings are another fit
  expect_false(identical(id_of(relabelled), id_of(loose)))
  expect_false(identical(id_of(stricter), id_of(loose)))
})


test_that("a partial refit with other settings keeps the recipe of the rows it left alone", {
  se <- .fit_hill(se_id, "mix", slicing_values = "GR", fit_args = list(pcutoff = 0.1))
  se <- .fit_hill(se, "mix", slicing_values = "RV", fit_args = list(pcutoff = 0.5))

  gr_id <- .hill_ids(se, "GR")
  rv_id <- .hill_ids(se, "RV")
  expect_length(gr_id, 1L)
  expect_length(rv_id, 1L)
  expect_false(identical(gr_id, rv_id))

  # both recipes are kept; before fit_id the second overwrote the first under "mix"
  expect_setequal(names(get_SE_provenance(se)$fits$hill), c(gr_id, rv_id))
  expect_identical(get_fit_recipe(se, "hill", gr_id)$fit_args, list(pcutoff = 0.1))
  expect_identical(get_fit_recipe(se, "hill", rv_id)$fit_args, list(pcutoff = 0.5))
  expect_identical(get_fit_provenance(se, "hill", fit_id = rv_id)$slicing_values, "RV")

  # asking by fit_source alone is ambiguous, and says so
  expect_error(get_fit_provenance(se, "hill", "mix"), "give the fit_id")

  summary <- summarize_fit_provenance(se)
  expect_setequal(summary$fit_id, c(gr_id, rv_id))
  expect_identical(unique(summary$fit_source), "mix")
  expect_true(all(c("fit_id", "solver") %in% names(summary)))

  checked <- check_reproducibility(se)
  expect_equal(NROW(checked), 2L)
  expect_true(all(checked$reproducible))
  expect_setequal(checked$fit_id, c(gr_id, rv_id))
})


test_that("refitting every slice with new settings drops the recipe nothing refers to", {
  se <- .fit_hill(se_id, "x", fit_args = list(pcutoff = 0.1))
  se <- .fit_hill(se, "x", fit_args = list(pcutoff = 0.5))

  expect_length(get_SE_provenance(se)$fits$hill, 1L)
  expect_identical(get_fit_provenance(se, "hill", "x")$fit_args, list(pcutoff = 0.5))
})


test_that("refit_from_provenance works per fit_id and leaves the other recipe alone", {
  se <- .fit_hill(se_id, "mix", slicing_values = "GR", fit_args = list(pcutoff = 0.1))
  se <- .fit_hill(se, "mix", slicing_values = "RV", fit_args = list(pcutoff = 0.5))
  gr_id <- .hill_ids(se, "GR")
  rv_id <- .hill_ids(se, "RV")

  refit <- suppressWarnings(refit_from_provenance(se, "hill", fit_id = gr_id))
  expect_identical(.hill_ids(refit, "GR"), gr_id)
  expect_identical(.hill_ids(refit, "RV"), rv_id)
  expect_setequal(names(get_SE_provenance(refit)$fits$hill), c(gr_id, rv_id))
  expect_equal(gDRutils:::.fit_rows(refit, "hill", "mix", gr_id),
               gDRutils:::.fit_rows(se, "hill", "mix", gr_id))
})


test_that("summary rows carry the fit_id and the summary recipe is keyed by it", {
  summary_fn <- .as_user(function(dt) list(n = NROW(dt)))
  se <- suppressWarnings(apply_fit(
    se_id, fit_drug_response_metrics, "single-agent",
    output_assay = "hill", fit_source = "s", summary_fn = summary_fn, summary_assay = "hill_n"
  ))
  id <- names(get_SE_provenance(se)$fits$hill)
  summary_ids <- unique(as.character(gDRutils:::.assay_dt(se, "hill_n")[["fit_id"]]))

  expect_identical(summary_ids, id)
  expect_identical(names(get_SE_provenance(se)$fits$hill_n), id)
  expect_identical(get_fit_provenance(se, "hill_n", "s")$fit_id, id)
})


test_that("objects from before fit_id read as if fit_source were the fit_id", {
  se_new <- .fit_hill(se_id, "loose", fit_args = list(pcutoff = 0.1))

  # the layout of schema 1: recipes keyed by fit_source, rows without a fit_id
  legacy <- se_new
  recipe <- get_SE_provenance(se_new)$fits$hill[[1]]
  recipe$fit_id <- NULL
  recipe$fit_source <- NULL
  recipe$solver <- NULL
  S4Vectors::metadata(legacy)$provenance <- list(
    schema_version = 1L, fits = list(hill = list(loose = recipe)), steps = list()
  )
  rows <- gDRutils:::.assay_dt(se_new, "hill")
  rows[, c("fit_id", "fit_message") := NULL]
  legacy <- persist_fit_assay(legacy, rows, "replace", "hill", "row", "column", "fit_source")

  expect_identical(get_fit_provenance(legacy, "hill", "loose")$fit_id, "loose")
  expect_identical(get_fit_recipe(legacy, "hill", "loose")$fit_source, "loose")
  expect_identical(summarize_fit_provenance(legacy)$fit_id, "loose")
  expect_equal(NROW(gDRutils:::.fit_rows(legacy, "hill", "loose", "loose")), NROW(rows))

  # a fit read back from schema 1 still reproduces, and is named by its hash afterwards
  checked <- suppressWarnings(check_reproducibility(legacy))
  expect_true(checked$reproducible)
  expect_identical(checked$fit_id, "loose")

  # recording a new fit names the old recipe and moves the object to schema 2
  next_se <- .fit_hill(legacy, "other", fit_args = list(pcutoff = 0.5))
  prov <- get_SE_provenance(next_se)
  expect_identical(prov$schema_version, 2L)
  expect_identical(prov$fits$hill$loose$fit_id, "loose")
  expect_identical(prov$fits$hill$loose$fit_source, "loose")
})


test_that("rows written without apply_fit stand for their recipe by fit_source", {
  se <- se_id
  rows <- data.table::data.table(
    row = rownames(se)[1], column = colnames(se)[1], fit_source = "PRISM",
    normalization_type = "RV", x_mean = 0.5
  )
  se <- persist_fit_assay(se, rows, "merge", "prism_like", "row", "column", "fit_source")
  se <- record_fit_provenance(se, "prism_like", "PRISM", function(se) se)
  expect_identical(get_fit_provenance(se, "prism_like", "PRISM")$fit_id, "PRISM")
  expect_identical(names(get_SE_provenance(se)$fits$prism_like), "PRISM")

  # and a later apply_fit into the same assay does not drop it
  se <- suppressWarnings(apply_fit(
    se, fit_drug_response_metrics, "single-agent", output_assay = "prism_like",
    fit_source = "gDR"
  ))
  expect_true("PRISM" %in% names(get_SE_provenance(se)$fits$prism_like))
  expect_equal(length(get_SE_provenance(se)$fits$prism_like), 2L)
})


test_that("the solver is recorded by its registry name", {
  se <- .fit_hill(se_id, "default")
  expect_identical(get_fit_provenance(se, "hill", "default")$solver, "drc")
  expect_identical(summarize_fit_provenance(se)$solver, "drc")

  explicit <- .fit_hill(se_id, "explicit", fit_args = list(solver = solver_drc))
  expect_identical(get_fit_provenance(explicit, "hill", "explicit")$solver, "drc")

  mine <- .as_user(function(conc, x, model, control = list()) {
    solver_drc(conc, x, model, control)
  })
  custom <- .fit_hill(se_id, "custom", fit_args = list(solver = mine))
  expect_identical(get_fit_provenance(custom, "hill", "custom")$solver, "custom")

  register_solver("test_fit_id_solver", mine, "a solver registered for a test")
  on.exit(rm(list = "test_fit_id_solver", envir = gDRutils:::.solver_registry), add = TRUE)
  registered <- .fit_hill(se_id, "registered", fit_args = list(solver = mine))
  expect_identical(get_fit_provenance(registered, "hill", "registered")$solver,
                   "test_fit_id_solver")

  # a solver changes the fit, so it changes the fit_id
  expect_false(identical(names(get_SE_provenance(custom)$fits$hill),
                         names(get_SE_provenance(se)$fits$hill)))
})


test_that("the model text no longer names a default solver", {
  expect_false(grepl("default solver_drc", get_fit_fn_info_model <- attr(
    fit_drug_response_metrics, "fit_fn_info"
  )$model))
})


test_that("fit_message holds the failure message of the solver, and is NA otherwise", {
  failing <- .as_user(function(conc, x, model, control = list()) {
    list(params = c(h = NA_real_, x_inf = NA_real_, x_0 = NA_real_, ec50 = NA_real_),
         rss = NA_real_, converged = FALSE, message = "the optimiser gave up")
  })
  dt <- data.table::data.table(
    Concentration = 10 ^ seq(-3, 1, length.out = 8),
    x = c(0.98, 0.95, 0.88, 0.70, 0.45, 0.25, 0.14, 0.10),
    x_std = 0.02,
    normalization_type = "RV"
  )
  bad <- fit_drug_response_metrics(dt, solver = failing)
  expect_identical(bad$fit_type, "DRCInvalidFitResult")
  expect_identical(bad$fit_message, "the optimiser gave up")

  throws <- .as_user(function(conc, x, model, control = list()) stop("it threw"))
  expect_identical(fit_drug_response_metrics(dt, solver = throws)$fit_message, "it threw")

  good <- fit_drug_response_metrics(dt)
  expect_true(is.na(good$fit_message))
  expect_true(is.character(good$fit_message))

  # constant fits and series with too few points have no solver message
  expect_true(is.na(fit_drug_response_metrics(dt[1:2])$fit_message))

  se <- .fit_hill(se_id, "failing", fit_args = list(solver = failing))
  rows <- gDRutils:::.assay_dt(se, "hill")
  expect_true(all(rows$fit_type[!is.na(rows$fit_type)] %in%
                    c("DRCInvalidFitResult", "DRCTooFewPointsToFit")))
  expect_true("the optimiser gave up" %in% rows$fit_message)
})


test_that("flatten leaves fit_id and fit_message out of the wide table, unless recipes differ", {
  tbl <- function(ids) {
    data.table::data.table(
      drug = "d", normalization_type = c("RV", "GR"), fit_source = "gDR",
      fit_type = "DRC4pHillFitModel", xc50 = c(1, 2),
      fit_id = ids, fit_message = c(NA_character_, "failed")
    )
  }
  groups <- c("normalization_type", "fit_source") # as convert_se_assay_to_dt() calls it
  wide <- c("fit_type", "xc50")

  one <- flatten(tbl(c("A", "A")), groups, wide)
  without <- flatten(tbl(c("A", "A"))[, c("fit_id", "fit_message") := NULL], groups, wide)
  expect_identical(one, without)
  expect_false(any(c("fit_id", "fit_message") %in% names(one)))

  several <- flatten(tbl(c("A", "B")), groups, wide)
  expect_true(all(c("RV_fit_id", "GR_fit_id") %in% names(several)))
  expect_identical(several$RV_fit_id, "A")
  expect_identical(several$GR_fit_id, "B")
  expect_identical(several$GR_fit_message, "failed")
})


test_that("objects fitted by different recipes, or by none, merge into the same treatments", {
  other_fn <- .as_user(function(dt) list(x_mean = mean(dt$x, na.rm = TRUE)))
  a <- .fit_hill(se_id, "gDR")
  b <- suppressWarnings(apply_fit(
    se_id, other_fn, "single-agent", output_assay = "hill", fit_source = "gDR"
  ))
  reference <- merge_SE(list(a = a, b = a))

  # a fit_id that differs between inputs is a value of the rows, not a treatment
  merged <- merge_SE(list(a = a, b = b))
  expect_identical(dim(merged), dim(reference))
  expect_true("fit_id" %in% names(convert_se_assay_to_dt(merged, "hill")))

  # and so is its absence in an input that was not fitted
  one_only <- merge_SE(list(fitted = a, bare = se_id))
  expect_identical(dim(one_only), dim(reference))
  expect_true(all(c("fit_id", "fit_message") %in% get_header("metrics_results")))
})
