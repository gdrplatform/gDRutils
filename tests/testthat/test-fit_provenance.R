suppressWarnings(
  se_prov <- qs2::qs_read(
    system.file("testdata/finalMAE_small.qs2", package = "gDRtestData")
  )[["single-agent"]]
)

# Tests run in an environment that inherits from the gDRutils namespace, so a
# function defined here would look like package code. Functions standing in for
# user code are therefore rebound to the global environment, optionally with
# values they capture, as a closure written in a session would be.
.as_user_fn <- function(fn, values = list()) {
  environment(fn) <- list2env(values, parent = globalenv())
  fn
}

mean_fn <- .as_user_fn(function(dt) list(x_mean = mean(dt$x, na.rm = TRUE)))

.rows <- function(se, assay, fit_source) .fit_rows(se, assay, fit_source)


test_that("apply_fit records a package function by reference with its fit_args", {
  se <- suppressWarnings(apply_fit(
    se_prov, fit_drug_response_metrics, "single-agent",
    output_assay = "hill", fit_source = "loose", fit_args = list(pcutoff = 0.1)
  ))
  recipe <- get_fit_provenance(se, "hill", "loose")

  expect_identical(recipe$called_by, "apply_fit")
  expect_identical(recipe$fn$name, "fit_drug_response_metrics")
  expect_identical(recipe$fn$package, "gDRutils")
  expect_identical(recipe$fn$version, as.character(utils::packageVersion("gDRutils")))
  expect_null(recipe$fn$source)
  expect_true("drc" %in% names(recipe$fn$dependencies))
  expect_match(recipe$fn$info$model, "LL.3u")
  expect_identical(recipe$fit_args, list(pcutoff = 0.1))
  expect_identical(recipe$input_assay, "Averaged")
  expect_identical(recipe$output_assay, "hill")
  expect_identical(recipe$slicing_values, c("GR", "RV"))
  expect_identical(recipe$run$n_rows, NROW(.rows(se, "hill", "loose")))
  expect_identical(get_SE_provenance(se)$schema_version, 1L)
})


test_that("fit_args reach the fit function and give the same result as a closure", {
  closure_fn <- .as_user_fn(
    function(dt) fit_drug_response_metrics(dt, pcutoff = p),
    values = list(p = 0.5)
  )
  se <- suppressWarnings(apply_fit(
    se_prov, fit_drug_response_metrics, "single-agent",
    output_assay = "hill", fit_source = "args", fit_args = list(pcutoff = 0.5)
  ))
  se <- suppressWarnings(apply_fit(
    se, closure_fn, "single-agent", output_assay = "hill", fit_source = "closure"
  ))

  by_args <- .rows(se, "hill", "args")
  by_closure <- .rows(se, "hill", "closure")
  by_args$fit_source <- by_closure$fit_source <- NULL
  expect_equal(by_args, by_closure)
})


test_that("a closure is stored as source with the values and helpers it uses", {
  user_fn <- .as_user_fn(
    function(dt) list(score = helper(dt$x) + offset),
    values = list(offset = 2, scale = 10)
  )
  helper <- function(x) mean(x, na.rm = TRUE) * scale
  environment(helper) <- environment(user_fn)
  assign("helper", helper, envir = environment(user_fn))

  se <- apply_fit(se_prov, user_fn, "single-agent", output_assay = "custom", fit_source = "user")
  fn <- get_fit_provenance(se, "custom", "user")$fn

  expect_null(fn$package)
  expect_match(fn$source, "helper(dt$x) + offset", fixed = TRUE)
  expect_identical(fn$captured, list(offset = 2, scale = 10))
  expect_match(fn$helpers$helper, "scale", fixed = TRUE)
  expect_null(fn$not_stored)
})


test_that("an inline function is recorded without a name", {
  se <- apply_fit(
    se_prov, function(dt) list(n = NROW(dt)), "single-agent",
    output_assay = "custom", fit_source = "inline"
  )
  fn <- get_fit_provenance(se, "custom", "inline")$fn
  expect_null(fn$name)
  expect_match(fn$source, "NROW(dt)", fixed = TRUE)
  expect_match(summarize_fit_provenance(se)$fn, "anonymous function (source stored)", fixed = TRUE)
})


test_that("a closure built inside a package function keeps the package and its settings", {
  factory_env <- new.env(parent = asNamespace("gDRutils"))
  factory_env$q <- 0.3
  pkg_closure <- eval(quote(function(dt) list(x_q = stats::quantile(dt$x, q, names = FALSE))),
                      envir = factory_env)

  se <- apply_fit(se_prov, pkg_closure, "single-agent", output_assay = "custom", fit_source = "q")
  fn <- get_fit_provenance(se, "custom", "q")$fn
  expect_identical(fn$package, "gDRutils")
  expect_identical(fn$captured, list(q = 0.3))
  expect_false(is.null(fn$source))

  res <- check_reproducibility(se, allow_source = TRUE)
  expect_true(res$reproducible)
})


test_that("values that cannot be stored are listed and reported", {
  big <- seq_len(1e5)
  user_fn <- .as_user_fn(function(dt) list(n = NROW(dt) + length(big)), values = list(big = big))
  expect_warning(
    se <- apply_fit(se_prov, user_fn, "single-agent", output_assay = "custom", fit_source = "big"),
    class = "gDRutils_provenance_warning"
  )
  expect_identical(get_fit_provenance(se, "custom", "big")$fn$not_stored, "big")

  ref_fn <- .as_user_fn(function(dt, ref) list(n = NROW(dt) + length(ref)))
  expect_warning(
    se <- apply_fit(se_prov, ref_fn, "single-agent", output_assay = "custom",
                    fit_source = "args", fit_args = list(ref = big)),
    class = "gDRutils_provenance_warning"
  )
  recipe <- get_fit_provenance(se, "custom", "args")
  expect_identical(recipe$args_not_stored, "ref")
  expect_null(recipe[["fit_args"]])
})


test_that("re-fitting replaces the recipe, and merge = 'replace' drops the others", {
  se <- apply_fit(se_prov, mean_fn, "single-agent", output_assay = "custom", fit_source = "a")
  se <- apply_fit(se, mean_fn, "single-agent", output_assay = "custom", fit_source = "b")
  se <- apply_fit(se, mean_fn, "single-agent", output_assay = "custom", fit_source = "a",
                  fit_args = list())
  expect_identical(names(get_SE_provenance(se)$fits$custom), c("a", "b"))

  se <- apply_fit(se, mean_fn, "single-agent", output_assay = "custom", fit_source = "c",
                  merge = "replace")
  expect_identical(names(get_SE_provenance(se)$fits$custom), "c")
})


test_that("a summary assay points at the recipe of the fit that wrote it", {
  summary_fn <- .as_user_fn(function(fit_dt) list(mean_all = mean(fit_dt$x_mean)))
  se <- apply_fit(
    se_prov, mean_fn, "single-agent", output_assay = "custom", fit_source = "s",
    summary_fn = summary_fn, summary_assay = "custom_summary"
  )
  prov <- get_SE_provenance(se)
  expect_identical(prov$fits$custom_summary$s, list(summary_of = "custom"))

  recipe <- get_fit_provenance(se, "custom_summary", "s")
  expect_identical(recipe$output_assay, "custom")
  expect_identical(recipe$summary_fn$name, "summary_fn")
  expect_true(check_reproducibility(se, allow_source = TRUE)$reproducible)
})


test_that("apply_fits records one recipe per output assay", {
  sd_fn <- .as_user_fn(function(dt) list(x_sd = stats::sd(dt$x) * k), values = list(k = 1))
  se <- apply_fits(
    se_prov,
    fit_fns = list(mean_assay = mean_fn, sd_assay = sd_fn),
    data_type = "single-agent",
    fit_source = "pair",
    fit_args = list(sd_assay = list())
  )
  prov <- get_SE_provenance(se)
  expect_setequal(names(prov$fits), c("mean_assay", "sd_assay"))
  expect_identical(prov$fits$sd_assay$pair$called_by, "apply_fits")
  expect_identical(prov$fits$sd_assay$pair$fn$name, "sd_fn")
  expect_identical(prov$fits$sd_assay$pair$fn$captured, list(k = 1))

  res <- check_reproducibility(se, allow_source = TRUE)
  expect_identical(NROW(res), 2L)
  expect_true(all(res$reproducible))
})


test_that("apply_fits records the function behind a shared pre-computation", {
  shared_fn <- .as_user_fn(function(dt) {
    m <- mean(dt$x, na.rm = TRUE)
    list(a_assay = list(a = m), b_assay = list(b = m * 2))
  })
  se <- apply_fits(
    se_prov,
    fit_fns = list(a_assay = shared_fn, b_assay = shared_fn),
    data_type = "single-agent",
    fit_source = "shared"
  )
  # one call fills both assays, so each slice gives exactly one row per assay
  expect_identical(NROW(.rows(se, "a_assay", "shared")), NROW(.rows(se, "b_assay", "shared")))
  expect_identical(NROW(.rows(se, "a_assay", "shared")), 2L * NROW(se_prov) * NCOL(se_prov))
  prov <- get_SE_provenance(se)
  expect_identical(prov$fits$a_assay$shared$fn$name, "shared_fn")
  expect_identical(prov$fits$b_assay$shared$fn$name, "shared_fn")
  expect_true(all(check_reproducibility(se, allow_source = TRUE)$reproducible))
})


test_that("check_reproducibility confirms, refuses and detects", {
  se <- suppressWarnings(apply_fit(
    se_prov, fit_drug_response_metrics, "single-agent",
    output_assay = "hill", fit_source = "pkg"
  ))
  se <- apply_fit(se, mean_fn, "single-agent", output_assay = "custom", fit_source = "user")

  res <- check_reproducibility(se)
  expect_true(res[fit_source == "pkg"]$reproducible)
  # code stored as source is not evaluated without consent
  expect_true(is.na(res[fit_source == "user"]$reproducible))
  expect_match(res[fit_source == "user"]$detail, "allow_source")

  res <- check_reproducibility(se, assay = "custom", allow_source = TRUE)
  expect_identical(res$fit_source, "user")
  expect_true(res$reproducible)

  # tamper with the stored results: the check must notice
  dt <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, "custom"), row.field = "row", column.field = "column"
  ))
  dt$x_mean <- dt$x_mean + 1
  se_bad <- persist_fit_assay(se, dt, "replace", "custom", "row", "column", "fit_source")
  res <- check_reproducibility(se_bad, assay = "custom", allow_source = TRUE)
  expect_false(res$reproducible)
})


test_that("refit_from_provenance warns when the recorded function changed", {
  se <- suppressWarnings(apply_fit(
    se_prov, fit_drug_response_metrics_4p, "single-agent",
    output_assay = "custom", fit_source = "v"
  ))
  prov <- get_SE_provenance(se)

  prov$fits$custom$v$fn$version <- "0.0.1"
  S4Vectors::metadata(se)$provenance <- prov
  expect_warning(
    suppressWarnings(refit_from_provenance(se, "custom", "v"), classes = "simpleWarning"),
    "recorded with version 0.0.1",
    class = "gDRutils_provenance_warning"
  )

  prov$fits$custom$v$fn$version <- as.character(utils::packageVersion("gDRutils"))
  prov$fits$custom$v$fn$hash <- "0000000000000000"
  S4Vectors::metadata(se)$provenance <- prov
  expect_warning(
    suppressWarnings(refit_from_provenance(se, "custom", "v"), classes = "simpleWarning"),
    "differs from the recorded one",
    class = "gDRutils_provenance_warning"
  )
})


test_that("the provenance is plain data, small, and survives qs2 and merge_SE", {
  se <- suppressWarnings(apply_fit(
    se_prov, fit_drug_response_metrics, "single-agent",
    output_assay = "Metrics", fit_source = "gDR"
  ))
  se <- apply_fit(se, mean_fn, "single-agent", output_assay = "custom", fit_source = "user")
  prov <- get_SE_provenance(se)

  expect_true(.is_plain_data(prov))
  # a recipe stays around a kilobyte once compressed
  expect_lt(length(qs2::qs_serialize(prov$fits$Metrics$gDR)), 3000)

  f <- tempfile(fileext = ".qs2")
  qs2::qs_save(se, f)
  expect_identical(get_SE_provenance(qs2::qs_read(f)), prov)

  merged <- merge_metadata(list(first = se, second = se_prov), "provenance")$provenance
  expect_identical(names(merged$sources), "first")
  expect_identical(merged$sources$first, prov)
  expect_identical(
    get_fit_provenance(se, "custom", "user"),
    .find_recipe(merged$sources$first, "custom", "user")
  )
})


test_that("merge_SE keeps the recipes of each input under its data_source", {
  se <- apply_fit(se_prov, mean_fn, "single-agent", output_assay = "custom", fit_source = "user")
  merged <- suppressWarnings(merge_SE(list(a = se, b = se)))

  expect_identical(names(get_SE_provenance(merged)$sources), c("a", "b"))
  expect_identical(get_fit_provenance(merged, "custom", "user", data_source = "b")$fn$name,
                   "mean_fn")
  summary <- summarize_fit_provenance(merged)
  expect_identical(summary$data_source, c("a", "b"))
  expect_error(get_fit_provenance(merged, "custom", "user", data_source = "c"), "data_source 'c'")
})


test_that("set_fit_fn_info attaches a description and a runtime hook", {
  described <- set_fit_fn_info(
    mean_fn,
    model = "mean response",
    description = "Mean of the averaged response.",
    reference = "none",
    env = function() list(tool = "1.0")
  )
  se <- apply_fit(se_prov, described, "single-agent", output_assay = "custom", fit_source = "d")
  fn <- get_fit_provenance(se, "custom", "d")$fn
  expect_identical(fn$info, list(model = "mean response",
                                 description = "Mean of the averaged response.",
                                 reference = "none"))
  expect_identical(fn$env, list(tool = "1.0"))
  # the description does not change the recorded code
  expect_identical(fn$hash, get_fit_provenance(
    apply_fit(se_prov, mean_fn, "single-agent", output_assay = "custom", fit_source = "d"),
    "custom", "d"
  )$fn$hash)

  expect_error(set_fit_fn_info(mean_fn, model = 1), "model")
  expect_error(set_fit_fn_info(mean_fn, env = function(x) x), "env")
})


test_that("getters handle objects without provenance and unknown fits", {
  expect_null(get_SE_provenance(se_prov))
  expect_identical(NROW(check_reproducibility(se_prov)), 0L)
  expect_identical(NROW(summarize_fit_provenance(se_prov)), 0L)
  se <- apply_fit(se_prov, mean_fn, "single-agent", output_assay = "custom", fit_source = "a")
  expect_error(get_fit_provenance(se, "custom", "nope"), "no recipe recorded")
})


test_that(".qualified_packages and .is_plain_data cover their edge cases", {
  expect_setequal(.qualified_packages(quote({
    data.table::setDT(x)
    gDRutils:::.drop_null(y)
    mean(z)
  })), c("data.table", "gDRutils"))
  expect_true(.is_plain_data(list(a = 1, b = list(c = "x"), d = data.frame(e = 1))))
  expect_false(.is_plain_data(list(f = function() 1)))
  expect_false(.is_plain_data(new.env()))
})
