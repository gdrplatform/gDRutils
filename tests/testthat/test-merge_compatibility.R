suppressWarnings(
  se_mc <- qs2::qs_read(
    system.file("testdata/finalMAE_small.qs2", package = "gDRtestData")
  )[["single-agent"]]
)

# the same data under another name for the drug column, with identifiers saying so
.with_drug_id <- function(se, name = "DrugID") {
  rd <- SummarizedExperiment::rowData(se)
  names(rd)[names(rd) == "Gnumber"] <- name
  SummarizedExperiment::rowData(se) <- rd
  ids <- get_SE_identifiers(se)
  ids$drug <- name
  suppressWarnings(set_SE_identifiers(se, ids))
}

.as_user_mc <- function(fn) {
  environment(fn) <- list2env(list(), parent = globalenv())
  fn
}

.fit_mc <- function(se, fit_source, fn = fit_drug_response_metrics, fit_args = list()) {
  suppressWarnings(apply_fit(
    se, fn, "single-agent", output_assay = "hill", fit_source = fit_source, fit_args = fit_args
  ))
}


test_that("merge_compatibility finds nothing to report for identical inputs", {
  report <- merge_compatibility(list(a = se_mc, b = se_mc))
  expect_s3_class(report, "merge_compatibility")
  expect_identical(NROW(report$issues), 0L)
  expect_setequal(report$identifiers$data_source, c("a", "b"))
  expect_false(any(report$identifiers$differs))
  expect_true(all(report$assays$present))
  expect_setequal(report$normalization_types$normalization_type, c("GR", "RV"))
  expect_output(print(report), "no differences")
})


test_that("merge_compatibility reports differing identifiers, assays and an input without provenance", {
  custom <- .with_drug_id(se_mc)
  fewer <- se_mc
  SummarizedExperiment::assays(fewer) <- SummarizedExperiment::assays(fewer)[
    setdiff(SummarizedExperiment::assayNames(fewer), "Controls")
  ]
  no_prov <- se_mc
  S4Vectors::metadata(no_prov)$provenance <- NULL

  report <- merge_compatibility(list(ref = se_mc, custom = custom, fewer = fewer, bare = no_prov))
  topics <- report$issues$topic
  expect_true(all(report$issues$severity == "difference"))
  expect_true("identifiers" %in% topics)
  expect_true("assays" %in% topics)
  expect_true("provenance" %in% topics)
  expect_identical(report$identifiers[report$identifiers$differs, ]$data_source, "custom")
  expect_identical(report$identifiers[report$identifiers$differs, ]$identifier, "drug")
  expect_false(report$assays[report$assays$data_source == "fewer" &
                               report$assays$assay == "Controls", ]$present)
  expect_output(print(report), "DrugID")
})


test_that("the same fit_source made by different methods is a contradiction; other settings are not", {
  other_fn <- .as_user_mc(function(dt) list(x_mean = mean(dt$x, na.rm = TRUE)))
  a <- .fit_mc(se_mc, "gDR")
  same_method <- .fit_mc(se_mc, "gDR", fit_args = list(pcutoff = 0.5))
  other_method <- .fit_mc(se_mc, "gDR", fn = other_fn)
  other_label <- .fit_mc(se_mc, "mine", fn = other_fn)

  settings <- merge_compatibility(list(a = a, b = same_method))
  expect_false("contradiction" %in% settings$issues$severity)
  expect_match(settings$issues$detail, "other settings")

  method <- merge_compatibility(list(a = a, b = other_method))
  expect_identical(method$issues[method$issues$severity == "contradiction", ]$topic, "fit_recipes")
  expect_match(method$issues[method$issues$severity == "contradiction", ]$detail, "fitted differently")

  expect_false("contradiction" %in% merge_compatibility(list(a = a, b = other_label))$issues$severity)
  expect_identical(NROW(merge_compatibility(list(a = a, b = a))$issues), 0L)
})


test_that("fits of different kinds under one label are a difference, not a contradiction", {
  # e.g. the Metrics of a combination experiment come from another step than those of a
  # single-agent experiment, both labelled "gDR"
  a <- .fit_mc(se_mc, "gDR")
  step_fn <- .as_user_mc(function(se) se)
  b <- record_fit_provenance(se_mc, "hill", "gDR", step_fn)

  report <- merge_compatibility(list(a = a, b = b))
  expect_false("contradiction" %in% report$issues$severity)
  expect_true(any(grepl("different steps", report$issues$detail)))
})


test_that("merge_SE stops on a contradiction unless allow_mismatch = TRUE", {
  other_fn <- .as_user_mc(function(dt) list(x_mean = mean(dt$x, na.rm = TRUE)))
  a <- .fit_mc(se_mc, "gDR")
  b <- .fit_mc(se_mc, "gDR", fn = other_fn)

  expect_error(merge_SE(list(a = a, b = b)), "contradict each other.*fitted differently")
  expect_error(merge_SE(list(a = a, b = b)), "allow_mismatch = TRUE")
  expect_warning(
    merged <- merge_SE(list(a = a, b = b), allow_mismatch = TRUE),
    "contradict each other"
  )
  report <- S4Vectors::metadata(merged)$merge_report
  expect_true(report$allow_mismatch)
  expect_true("contradiction" %in% report$issues$severity)

  # and the same fits merge without a word
  expect_silent(merge_SE(list(a = a, b = a)))
})


test_that("merge_SE keeps its merge_report in the metadata, as plain data", {
  merged <- merge_SE(list(a = se_mc, b = se_mc))
  report <- S4Vectors::metadata(merged)$merge_report

  expect_true(all(c("identifiers", "assays", "fit_recipes", "issues", "target_identifiers",
                    "allow_mismatch") %in% names(report)))
  expect_false(report$allow_mismatch)
  expect_identical(report$target_identifiers, get_SE_identifiers(se_mc))
  expect_true(all(vapply(report[c("identifiers", "issues")], is.data.frame, logical(1))))

  f <- tempfile(fileext = ".qs2")
  qs2::qs_save(merged, f)
  expect_identical(S4Vectors::metadata(qs2::qs_read(f))$merge_report, report)
})


test_that("merge_SE renames the identifier columns of the other inputs and leaves the session alone", {
  reset_env_identifiers()
  on.exit(reset_env_identifiers(), add = TRUE)
  custom <- .with_drug_id(se_mc)

  # the first input names the columns of the result
  merged_custom_first <- merge_SE(list(custom = custom, ref = se_mc))
  expect_true("DrugID" %in% names(SummarizedExperiment::rowData(merged_custom_first)))
  expect_false("Gnumber" %in% names(SummarizedExperiment::rowData(merged_custom_first)))
  expect_identical(get_SE_identifiers(merged_custom_first, "drug"), "DrugID")
  # the session identifiers are as they were
  expect_identical(get_env_identifiers("drug"), "Gnumber")

  merged_ref_first <- merge_SE(list(ref = se_mc, custom = custom))
  expect_true("Gnumber" %in% names(SummarizedExperiment::rowData(merged_ref_first)))
  expect_false("DrugID" %in% names(SummarizedExperiment::rowData(merged_ref_first)))
  expect_identical(get_env_identifiers("drug"), "Gnumber")

  # the same data under another name merges into the same treatments, not into twice as many
  baseline <- merge_SE(list(a = se_mc, b = se_mc))
  expect_identical(dim(merged_ref_first), dim(baseline))
  expect_identical(NROW(SummarizedExperiment::rowData(merged_custom_first)),
                   NROW(SummarizedExperiment::rowData(baseline)))
  merged_dt <- convert_se_assay_to_dt(merged_ref_first, "Averaged")
  expect_false(any(is.na(merged_dt[["Gnumber"]])))

  report <- S4Vectors::metadata(merged_ref_first)$merge_report
  expect_true("identifiers" %in% report$issues$topic)
})


test_that("merge_SE keeps every input under provenance$sources, with its metadata", {
  fitted <- .fit_mc(se_mc, "mine")
  bare <- se_mc
  S4Vectors::metadata(bare)$provenance <- NULL

  merged <- suppressWarnings(merge_SE(list(fitted = fitted, bare = bare)))
  sources <- get_SE_provenance(merged)$sources

  expect_identical(names(sources), c("fitted", "bare"))
  expect_false(isTRUE(sources$fitted$recorded))
  expect_false(sources$bare$recorded)
  expect_identical(sources$bare$fits, list())
  expect_identical(sources$fitted$metadata$identifiers, S4Vectors::metadata(fitted)$identifiers)
  expect_identical(sources$fitted$metadata$fit_parameters, get_SE_fit_parameters(fitted))
  expect_identical(sources$fitted$metadata$versions, get_SE_processing_metadata(fitted)$versions)
  expect_identical(get_fit_provenance(merged, "hill", "mine", data_source = "fitted")$fit_source, "mine")
})


test_that("merge_MAE passes additional_col_name and discard_keys on to merge_SE", {
  mae <- get_synthetic_data("finalMAE_small")
  exp_name <- names(mae)[1]

  by_se <- merge_SE(list(a = mae[[exp_name]], b = mae[[exp_name]]), additional_col_name = "origin")
  by_mae <- suppressWarnings(merge_MAE(list(a = mae, b = mae), additional_col_name = "origin"))
  by_default <- suppressWarnings(merge_MAE(list(a = mae, b = mae)))

  cols <- function(se) names(convert_se_assay_to_dt(se, "Metrics"))
  expect_identical(cols(by_mae[[exp_name]]), cols(by_se))
  expect_true("origin" %in% cols(by_mae[[exp_name]]))
  expect_false("origin" %in% cols(by_default[[exp_name]]))
  expect_true("data_source" %in% cols(by_default[[exp_name]]))

  # merge_SE checks both arguments, so a bad value only errors if merge_MAE passes it on
  expect_error(suppressWarnings(merge_MAE(list(a = mae, b = mae), discard_keys = 1)), "character")
  expect_error(suppressWarnings(merge_MAE(list(a = mae, b = mae), additional_col_name = 1)), "string")
})
