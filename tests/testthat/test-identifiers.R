test_that("get_env_identifiers and set_env_identifier work", {
  reset_env_identifiers()

  expect_error(get_env_identifiers("BOGUS", simplify = TRUE))
  expect_equal(get_env_identifiers("duration", simplify = TRUE), "Duration")
  expect_equal(get_env_identifiers(c("duration", "barcode"), simplify = FALSE),
               list(duration = "Duration", barcode = c("Barcode", "Plate")))

  set_env_identifier("cellline", "my_personal_cell_line_identifiers")
  expect_equal(get_env_identifiers("cellline", simplify = TRUE), "my_personal_cell_line_identifiers")
  vals <- IDENTIFIERS_LIST
  vals[["cellline"]] <- "my_personal_cell_line_identifiers"

  expect_equal(get_env_identifiers(simplify = TRUE), vals)
})


test_that("reset_env_identifiers works", {
  reset_env_identifiers()

  d <- get_env_identifiers(simplify = TRUE)
  set_env_identifier("duration", "TEST_DURATION")
  expect_equal(get_env_identifiers("duration", simplify = TRUE), "TEST_DURATION")

  ilist <- reset_env_identifiers()
  expect_equal(ilist, NULL)
  expect_equal(get_env_identifiers("duration", simplify = TRUE), "Duration")
})


test_that("get_prettified_identifiers works as expected", {
  obs <- get_prettified_identifiers(c("drug_name", "cellline_name"), simplify = FALSE)
  expect_equal(obs, c("Drug Name", "Cell Line Name"))

  obs <- get_prettified_identifiers(simplify = TRUE)
  expect_true(is(obs, "list"))
  expect_true(length(obs) > 1L)
})

test_that("get_env_identifier works with untreated_tag", {
  obs <- get_env_identifiers("untreated_tag", simplify = TRUE)
  expect_length(obs, 2)
  expect_equal(obs, c("vehicle", "untreated"))
})

test_that("get_SE_identifier works with untreated_tag", {
  se <- SummarizedExperiment()
  expect_warning(
    obs <- get_SE_identifiers(se, "untreated_tag", simplify = TRUE),
    "'se' was passed, but identifier 'untreated_tag' not found on se's identifiers"
  )
  expect_length(obs, 2)
  expect_equal(obs, c("vehicle", "untreated"))
})

test_that("get_expect_one_identifiers works as expected", {
  expect_true(all(get_expect_one_identifiers() %in% names(get_env_identifiers())))
  expect_true(length(get_expect_one_identifiers()) > 1L)
})

test_that("get_required_identifiers works as expected", {
  expect_true(all(get_required_identifiers() %in% names(get_env_identifiers())))
  expect_true(length(get_required_identifiers()) > 1L)
  expect_true(all(get_required_identifiers() %in% get_expect_one_identifiers()))
})

test_that("get_identifiers_dt works as expected", {
  dt <- yaml::read_yaml(system.file(package = "gDRutils", "identifier_descriptions.yaml"))
  expect_equal(get_identifiers_dt(), dt)
  expect_equal(get_identifiers_dt(k = "drug_name"), dt[["drug_name"]])
  expect_equal(get_identifiers_dt(get_description = TRUE),
               lapply(dt, function(i) {
                 i[["description"]]
                 }))
  expect_equal(get_identifiers_dt(get_example = TRUE),
               lapply(dt, function(i) {
                 i[["example"]]
                 }))
  expect_equal(get_identifiers_dt(k = "drug_name", get_description = TRUE),
               dt$drug_name$description)
  expect_equal(get_identifiers_dt(k = "drug_name", get_example = TRUE),
               dt$drug_name$example)
  expect_error(get_identifiers_dt(k = "some_drugs"),
               "Assertion on 'k %in% names(dt)' failed: Must be TRUE.",
               fixed = TRUE)
  expect_error(get_identifiers_dt(k = 2),
               "Assertion on 'k' failed: Must be of type 'string' (or 'NULL'), not 'double'.",
               fixed = TRUE)
})

test_that("'get_synonyms works", {
  expect_true(is.list(get_idfs_synonyms()))
  expect_true(length(get_idfs_synonyms()) > 0)
})

test_that("'update_synonyms works", {
  mdict <- list(duration = "time")
  dval <- gDRutils::get_env_identifiers("duration")
  iv <- c("Time", dval, "time")
  out <- update_idfs_synonyms(iv, dict = mdict)
  expect_identical(unique(out), dval)

  out <- update_idfs_synonyms(list(a = iv, b = iv), dict = mdict)
  expect_true(length(out) == 2)
  expect_identical(unique(out[[1]]), dval)
})

test_that("update_env_identifers_from_mae works as expected", {
  mae_idfs <- list()
  mae_idfs[["first_SE"]] <- gDRutils::get_env_identifiers()
  mae_idfs[["first_SE"]][["barcode"]] <- "my_barcode"

  mae_idfs1 <- mae_idfs
  mae_idfs1[["first_SE"]][["barcode"]] <- c("my_barcode", "my_plate")

  mae_idfs2 <- mae_idfs
  mae_idfs2[["second_SE"]] <- mae_idfs[["first_SE"]]

  mae_idfs3 <- mae_idfs1
  mae_idfs3[["second_SE"]] <- mae_idfs1[["first_SE"]]

  mae_idfs4 <- mae_idfs
  mae_idfs4[["second_SE"]] <- mae_idfs1[["first_SE"]]

  gDRutils::reset_env_identifiers()
  gDRutils::update_env_idfs_from_mae(mae_idfs)
  expect_equal(gDRutils::get_env_identifiers("barcode"), "my_barcode")

  gDRutils::reset_env_identifiers()
  gDRutils::update_env_idfs_from_mae(mae_idfs1)
  expect_equal(gDRutils::get_env_identifiers("barcode"), c("my_barcode", "my_plate"))

  gDRutils::reset_env_identifiers()
  gDRutils::update_env_idfs_from_mae(mae_idfs2)
  expect_equal(gDRutils::get_env_identifiers("barcode"), "my_barcode")

  gDRutils::reset_env_identifiers()
  gDRutils::update_env_idfs_from_mae(mae_idfs3)
  expect_equal(gDRutils::get_env_identifiers("barcode"), c("my_barcode", "my_plate"))

  gDRutils::reset_env_identifiers()
  expect_error(gDRutils::update_env_idfs_from_mae(mae_idfs4),
               "identical(mae_idfs[[1]], mae_idfs[[2]]) is not TRUE",
               fixed = TRUE)
  expect_error(gDRutils::update_env_idfs_from_mae("test"),
               "Assertion on 'mae_idfs' failed: Must be of type 'list', not 'character'.",
               fixed = TRUE)

})
