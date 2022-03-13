library(testthat); library(gDRutils)

test_that(".check_required_identifiers works as expected", {
  # Set up.
  nrow <- 5
  req_ids <- get_required_identifiers()
  nids <- length(req_ids)

  ids <- as.list(LETTERS[seq(nids)])
  names(ids) <- req_ids

  df <- as.data.frame(matrix(rep(0, nids * nrow), nrow = nrow, ncol = nids))
  names(df) <- ids

  # Single mapping all exist.
  obs <- gDRutils:::.check_required_identifiers(df, req_ids = req_ids, id_map = ids)
  expect_equal(obs, NULL)

  # Single mapping more ids than required exist.
  excess_map_ids <- ids
  excess_map_ids <- c(excess_map_ids, list(test = "excess"))
  obs <- gDRutils:::.check_required_identifiers(df, req_ids = req_ids, id_map = excess_map_ids)
  expect_equal(obs, NULL)

  # Single mapping more columns than required exist.
  fewer_map_ids <- ids
  fewer_req_ids <- req_ids[seq(length(req_ids) - 1)]
  fewer_map_ids <- fewer_map_ids[req_ids]
  obs <- gDRutils:::.check_required_identifiers(df, req_ids = fewer_req_ids, id_map = fewer_map_ids)
  expect_equal(obs, NULL)

  # Single mapping required id value does not exist.
  single_map_ids <- ids
  single_map_ids[[1]] <- "Cinderella"
  obs <- gDRutils:::.check_required_identifiers(df, req_ids = req_ids, id_map = single_map_ids)
  exp <- "specified value identifier(s): 'Cinderella' do not exist for standardized identifier(s): 'duration'"
  expect_equal(obs, exp)

  # Missing identifier.
  missing_map_ids <- ids
  missing_map_ids[[1]] <- NULL
  expect_error(gDRutils:::.check_required_identifiers(df, req_ids = req_ids, id_map = missing_map_ids),
    regex = sprintf("required identifiers: '%s' missing in 'id_map'", names(ids)[[1]]))

  # Polymapping. 
  poly_map_ids <- ids
  poly_map_ids[[1]] <- c("Cinderella", "Mulan")
  expect_error(gDRutils:::.check_required_identifiers(df, req_ids = req_ids, id_map = poly_map_ids),
    regex = sprintf("more than one identifier value found for the required identifiers: '%s'", names(poly_map_ids)[1]))
})

test_that(".check_polymapped_identifiers works as expected", {
  # Set up.
  nrow <- 5
  req_ids <- get_required_identifiers()
  nids <- length(req_ids)

  ids <- as.list(LETTERS[seq(nids)])
  names(ids) <- req_ids

  df <- as.data.frame(matrix(rep(0, nids * nrow), nrow = nrow, ncol = nids))
  names(df) <- ids

  # TODO
})

test_that("gDRutils:::.singlefy_polymapped_identifiers works as expected", {
  # Set up.

  gDRutils:::.singlefy_polymapped_identifiers()

  gDRutils:::.singlefy_polymapped_identifiers()
})

test_that("validate_identifiers works as expected", {
  # Set up.
  nrow <- 5
  req_ids <- get_required_identifiers()
  nids <- length(req_ids)

  ids <- as.list(LETTERS[seq(nids)])
  names(ids) <- req_ids

  df <- as.data.frame(matrix(rep(0, nids * nrow), nrow = nrow, ncol = nids))
  names(df) <- ids

  # Single mapping but does not exist.
  single_map_ids <- ids
  single_map_ids[1] <- "Cinderella"
  expect_error(validate_identifiers(df, identifiers = single_map_ids, req_ids = req_ids),
    regex = "")

  # TODO
  # Polymapping not
  expect_error(validate_identifiers(df, identifiers = ))
})
