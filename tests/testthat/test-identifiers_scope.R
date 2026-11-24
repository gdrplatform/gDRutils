test_that("local_identifiers restores the session on return and on error", {
  reset_env_identifiers()
  f <- function() {
    local_identifiers(list(drug = "DrugID", cellline = "Clid"))
    c(get_env_identifiers("drug"), get_env_identifiers("cellline"))
  }
  expect_identical(f(), c("DrugID", "Clid"))
  expect_identical(get_env_identifiers(), IDENTIFIERS_LIST)

  g <- function() {
    local_identifiers(list(drug = "DrugID"))
    stop("boom")
  }
  expect_error(g(), "boom")
  expect_identical(get_env_identifiers("drug"), "Gnumber")
})


test_that("scopes nest and keep the keys they do not set", {
  reset_env_identifiers()
  set_env_identifier("duration", "Time")
  on.exit(reset_env_identifiers(), add = TRUE)
  inner <- function() {
    local_identifiers(list(drug = "inner"))
    get_env_identifiers("drug")
  }
  outer <- function() {
    local_identifiers(list(drug = "outer"))
    c(inner(), get_env_identifiers("drug"), get_env_identifiers("duration"))
  }
  expect_identical(outer(), c("inner", "outer", "Time"))
  expect_identical(get_env_identifiers("drug"), "Gnumber")
  expect_identical(get_env_identifiers("duration"), "Time")
})


test_that("with_identifiers and with_SE_identifiers scope one expression", {
  reset_env_identifiers()
  expect_identical(with_identifiers(list(cellline = "Clid"), get_env_identifiers("cellline")), "Clid")
  se <- SummarizedExperiment::SummarizedExperiment()
  se <- set_SE_identifiers(se, list(drug = "DrugID"))
  expect_identical(with_SE_identifiers(se, get_env_identifiers("drug")), "DrugID")
  expect_identical(get_env_identifiers("drug"), "Gnumber")
  expect_error(local_identifiers(list(bogus = "x")), "bogus")
})


test_that("an SE's own identifiers are used to read it, and the session is left alone", {
  reset_env_identifiers()
  se <- get_synthetic_data("finalMAE_small")[["single-agent"]]
  ref <- convert_se_assay_to_dt(se, "Averaged")

  # the same data under other column names, with identifiers saying so
  custom <- se
  rd <- SummarizedExperiment::rowData(custom)
  names(rd)[names(rd) == "Gnumber"] <- "DrugID"
  SummarizedExperiment::rowData(custom) <- rd
  ids <- get_SE_identifiers(custom)
  ids$drug <- "DrugID"
  custom <- suppressWarnings(set_SE_identifiers(custom, ids))  # replaces the stored ones

  got <- convert_se_assay_to_dt(custom, "Averaged")
  expect_true("DrugID" %in% names(got))
  data.table::setnames(got, "DrugID", "Gnumber")
  expect_equal(got[, names(ref), with = FALSE], ref)
  expect_identical(get_env_identifiers("drug"), "Gnumber")

  # everything beneath the entry point sees the object's identifiers
  seen <- NULL
  local_mocked_bindings(.convert_se_assay_to_dt = function(...) {
    seen <<- get_env_identifiers("drug")
    stop("seen")
  })
  expect_error(convert_se_assay_to_dt(custom, "Averaged"), "seen")
  expect_identical(seen, "DrugID")
  expect_identical(get_env_identifiers("drug"), "Gnumber")
})


test_that("canonical_identifiers completes, resolves and marks assumed keys", {
  ids <- canonical_identifiers(list(barcode = c("Barcode", "Plate"), drug = "DrugID"),
                               data = c("Plate", "DrugID"))
  expect_identical(ids$barcode, "Plate")
  expect_identical(ids$drug, "DrugID")
  expect_identical(names(ids), names(IDENTIFIERS_LIST))
  expect_false("drug" %in% attr(ids, "assumed"))
  expect_true("cellline" %in% attr(ids, "assumed"))
  # with no alternative present, the list is kept rather than lost
  expect_identical(canonical_identifiers(list(barcode = c("Barcode", "Plate")), data = "x")$barcode,
                   c("Barcode", "Plate"))
  # a key that is no longer an identifier is left out and named
  ids <- canonical_identifiers(list(bogus = "x", drug = "DrugID"))
  expect_false("bogus" %in% names(ids))
  expect_identical(attr(ids, "unknown"), "bogus")
})


test_that("same_identifiers treats a resolved alternative as the same column", {
  expect_true(same_identifiers(list(barcode = "Barcode"), IDENTIFIERS_LIST))
  expect_true(same_identifiers(list(), IDENTIFIERS_LIST))
  res <- same_identifiers(list(drug = "DrugID", barcode = "Other"), IDENTIFIERS_LIST)
  expect_false(res)
  expect_identical(sort(attr(res, "differ")), c("barcode", "drug"))
})


test_that("a key the object does not store means the session's value, for comparing too", {
  reset_env_identifiers()
  on.exit(reset_env_identifiers(), add = TRUE)
  set_env_identifier("drug", "d1_id")
  a <- IDENTIFIERS_LIST[setdiff(names(IDENTIFIERS_LIST), "drug")]
  b <- IDENTIFIERS_LIST
  se_a <- set_SE_identifiers(SummarizedExperiment::SummarizedExperiment(), a)
  se_b <- set_SE_identifiers(SummarizedExperiment::SummarizedExperiment(), b)
  # they are read with different columns, so they must not compare equal
  expect_identical(with_SE_identifiers(se_a, get_env_identifiers("drug")), "d1_id")
  expect_identical(with_SE_identifiers(se_b, get_env_identifiers("drug")), "Gnumber")
  expect_false(same_identifiers(a, b))
  expect_error(update_env_idfs_from_mae(list(a, b)), "different identifiers: drug")
})


test_that("tags and alternatives compare as sets; well positions keep their order", {
  expect_true(same_identifiers(list(untreated_tag = c("untreated", "vehicle")), IDENTIFIERS_LIST))
  expect_true(same_identifiers(list(barcode = c("Plate", "Barcode")), IDENTIFIERS_LIST))
  expect_false(same_identifiers(list(well_position = c("WellColumn", "WellRow")), IDENTIFIERS_LIST))
})


test_that("update_env_idfs_from_mae warns about keys that are no longer identifiers", {
  reset_env_identifiers()
  on.exit(reset_env_identifiers(), add = TRUE)
  old <- c(IDENTIFIERS_LIST, list(retired_key = "x"))
  expect_warning(update_env_idfs_from_mae(list(old)), "retired_key")
  expect_warning(update_env_idfs_from_mae(list(old, old)), "retired_key")
  expect_identical(get_env_identifiers(), IDENTIFIERS_LIST)
})


test_that("identifier values must be character; a NULL stored on an SE counts as not stored", {
  reset_env_identifiers()
  expect_error(local_identifiers(list(drug2 = NULL)), "drug2|character")
  se <- set_SE_identifiers(SummarizedExperiment::SummarizedExperiment(),
                           list(drug2 = NULL, drug = "DrugID"))
  expect_identical(with_SE_identifiers(se, get_env_identifiers("drug2")), "Gnumber_2")
  expect_identical(with_SE_identifiers(se, get_env_identifiers("drug")), "DrugID")
})


test_that("local_identifiers at the top level warns that nothing restores it", {
  reset_env_identifiers()
  on.exit(reset_env_identifiers(), add = TRUE)
  expect_warning(local_identifiers(list(drug = "X"), .local_envir = globalenv()), "top level")
  expect_identical(get_env_identifiers("drug"), "X")
})


test_that("the keys convert_se_assay_to_dt uses come from the SE, not the session", {
  reset_env_identifiers()
  se <- get_synthetic_data("finalMAE_small")[["single-agent"]]
  ref <- convert_se_assay_to_dt(se, "Normalized")
  read_in_other_session <- function() {
    local_identifiers(list(concentration = "conc_elsewhere", masked_tag = "mask_elsewhere"))
    convert_se_assay_to_dt(se, "Normalized")
  }
  expect_equal(read_in_other_session(), ref)
})


test_that("validate_SE stops when the SE's identifiers name columns it does not have", {
  reset_env_identifiers()
  se <- get_synthetic_data("finalMAE_small")[["single-agent"]]
  expect_no_error(validate_SE(se))
  ids <- get_SE_identifiers(se)
  ids$drug <- "no_such_column"
  bad <- suppressWarnings(set_SE_identifiers(se, ids))
  expect_error(validate_SE(bad), "drug = 'no_such_column'")
})
