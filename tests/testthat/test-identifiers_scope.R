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
  expect_error(canonical_identifiers(list(bogus = "x")), "bogus")
})


test_that("same_identifiers treats a resolved alternative as the same column", {
  expect_true(same_identifiers(list(barcode = "Barcode"), IDENTIFIERS_LIST))
  expect_true(same_identifiers(list(), IDENTIFIERS_LIST))
  res <- same_identifiers(list(drug = "DrugID", barcode = "Other"), IDENTIFIERS_LIST)
  expect_false(res)
  expect_identical(sort(attr(res, "differ")), c("barcode", "drug"))
})
