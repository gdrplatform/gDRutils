test_that("the built-ins are registered and resolvable", {
  tbl <- get_fit_fns()
  expect_s3_class(tbl, "data.frame")
  expect_true(all(c("fit_drug_response_metrics",
                    "fit_drug_response_metrics_4p") %in% tbl$name))
  expect_identical(get_fit_fn("fit_drug_response_metrics"),
                   fit_drug_response_metrics)
})


test_that("data_type filters, and an unmatched one gives zero rows", {
  sa <- get_fit_fns(data_type = "single-agent")
  expect_true(NROW(sa) >= 2L)
  expect_true(all(grepl("single-agent", sa$data_type, fixed = TRUE)))

  none <- get_fit_fns(data_type = "no-such-profile")
  expect_identical(NROW(none), 0L)
  # the empty case must keep the columns, or a caller that filters then selects
  # fails only when nothing matches, which is the hardest time to notice
  expect_identical(names(none), names(sa))
})


test_that("a function may be registered for several data types", {
  withr::defer(rm("multi_type_fn", envir = gDRutils:::.fit_fn_env))
  register_fit_fn(
    "multi_type_fn",
    package = "gDRutils",
    fn_name = "fit_drug_response_metrics",
    data_type = c("single-agent", "combination")
  )
  expect_true("multi_type_fn" %in% get_fit_fns("single-agent")$name)
  expect_true("multi_type_fn" %in% get_fit_fns("combination")$name)
})


test_that("a function carrying no fit_fn_info lists as NA, not as empty text", {
  # NA rather than "" so a caller can tell "nobody described this" from
  # "described as nothing".
  row <- get_fit_fns()[get_fit_fns()$name == "fit_drug_response_metrics", ]
  expect_true(is.na(row$model))
  expect_true(is.na(row$description))
  expect_true(is.na(row$reference))
})


test_that("fit_fn_info is read into the table when a function carries it", {
  # Tested on the field reader rather than end to end: the registry resolves
  # through getExportedValue(), so an end-to-end case needs an *exported*
  # function carrying the attribute, and gDRutils has none. gDRcore's
  # musyc_summary_fn does, which is where the round trip is exercised.
  info <- attr(
    set_fit_fn_info(
      function(dt) list(x_mean = mean(dt$x, na.rm = TRUE)),
      model = "mean response",
      description = "Mean of the averaged response.",
      reference = "none"
    ),
    "fit_fn_info"
  )
  expect_identical(gDRutils:::.info_field(info, "model"), "mean response")
  expect_identical(gDRutils:::.info_field(info, "reference"), "none")
  expect_true(is.na(gDRutils:::.info_field(info, "absent_key")))
})


test_that("an unavailable predicate hides the function without an error", {
  withr::defer(rm("needs_python_fn", envir = gDRutils:::.fit_fn_env))
  register_fit_fn(
    "needs_python_fn",
    package = "gDRutils",
    fn_name = "fit_drug_response_metrics",
    data_type = "combination",
    available = function() FALSE
  )
  row <- get_fit_fns()[get_fit_fns()$name == "needs_python_fn", ]
  expect_false(row$available)
})


test_that("a predicate that errors counts as unavailable, not as a failure", {
  withr::defer(rm("broken_predicate_fn", envir = gDRutils:::.fit_fn_env))
  register_fit_fn(
    "broken_predicate_fn",
    package = "gDRutils",
    fn_name = "fit_drug_response_metrics",
    data_type = "combination",
    available = function() stop("no python here")
  )
  expect_silent(tbl <- get_fit_fns())
  expect_false(tbl[tbl$name == "broken_predicate_fn", ]$available)
})


test_that("an entry from a package that is not installed lists as unavailable", {
  withr::defer(rm("absent_pkg_fn", envir = gDRutils:::.fit_fn_env))
  register_fit_fn(
    "absent_pkg_fn",
    package = "notAnInstalledPackage",
    data_type = "single-agent"
  )
  # listing must not fail because one provider is missing
  expect_silent(tbl <- get_fit_fns())
  expect_false(tbl[tbl$name == "absent_pkg_fn", ]$available)
  # but asking for it directly is an error, not a silent NULL
  expect_error(get_fit_fn("absent_pkg_fn"), "does not export it")
})


test_that("registering the same name twice overwrites", {
  withr::defer(rm("dup_fn", envir = gDRutils:::.fit_fn_env))
  register_fit_fn("dup_fn", package = "gDRutils",
                  fn_name = "fit_drug_response_metrics",
                  data_type = "single-agent")
  register_fit_fn("dup_fn", package = "gDRutils",
                  fn_name = "fit_drug_response_metrics",
                  data_type = "combination")
  expect_identical(NROW(get_fit_fns()[get_fit_fns()$name == "dup_fn", ]), 1L)
  expect_true("dup_fn" %in% get_fit_fns("combination")$name)
  expect_false("dup_fn" %in% get_fit_fns("single-agent")$name)
})


test_that("an unknown name is an error that lists what is registered", {
  expect_error(get_fit_fn("no_such_fit_fn"), "Unknown fit function")
  expect_error(get_fit_fn("no_such_fit_fn"), "fit_drug_response_metrics")
})


test_that("arguments are validated", {
  expect_error(register_fit_fn(1, package = "gDRutils",
                               data_type = "single-agent"))
  expect_error(register_fit_fn("x", package = "gDRutils",
                               data_type = character(0)))
  expect_error(register_fit_fn("x", package = "gDRutils",
                               data_type = "single-agent",
                               available = function(a) TRUE))
  expect_error(get_fit_fns(data_type = 1))
})
