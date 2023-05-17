library(testthat)
library(gDRutils)

test_that(".check_required_identifiers works as expected", {
  # Set up.
  nrow <- 5
  req_ids <- get_required_identifiers()
  nids <- length(req_ids)

  ids <- as.list(LETTERS[seq(nids)])
  names(ids) <- req_ids

  df <- data.table::data.table(matrix(rep(0, nids * nrow), nrow = nrow, ncol = nids))
  colnames(df) <- unlist(ids)

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
  exp <- "specified value identifier(s): 'Cinderella' do not exist for standardized identifier(s): 'duration'\n"
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
    regex = sprintf("more than one identifier value found for required identifiers: '%s'", names(poly_map_ids)[1]))
})

test_that(".check_polymapped_identifiers works as expected", {
  # Set up.
  nrow <- 5
  exp_one_ids <- get_expect_one_identifiers()
  nids <- length(exp_one_ids)

  id_map <- as.list(LETTERS[seq(nids)])
  names(id_map) <- exp_one_ids

  df <- data.table::data.table(matrix(rep(0, nids * nrow), nrow = nrow, ncol = nids))
  colnames(df) <- unlist(id_map)

  # All singletons.
  obs <- gDRutils:::.check_polymapped_identifiers(df, exp_one_ids, id_map)
  expect_equal(obs, NULL)

  # Some polymappings.
  some_poly_map <- id_map
  some_poly_map[[1]] <- c(some_poly_map[[1]], "extra_item") 
  obs <- gDRutils:::.check_polymapped_identifiers(df, exp_one_ids, id_map = some_poly_map)
  expect_equal(obs, "more than one mapping for identifier(s): 'duration'\n")

  # All polymappings.
  all_poly_map <- id_map
  all_poly_map <- lapply(seq_along(all_poly_map), function(x) c(all_poly_map[[x]], "extra_item"))
  names(all_poly_map) <- names(id_map)
  obs <- gDRutils:::.check_polymapped_identifiers(df, exp_one_ids, id_map = all_poly_map)
  exp <- paste0(names(id_map), collapse = ", ")
  expect_equal(obs, sprintf("more than one mapping for identifier(s): '%s'\n", exp))
})

test_that("gDRutils:::.modify_polymapped_identifiers works as expected", {
  # Set up.
  nrow <- 5
  exp_one_ids <- get_expect_one_identifiers()
  nids <- length(exp_one_ids)

  id_map <- as.list(LETTERS[seq(nids)])
  names(id_map) <- exp_one_ids

  df <- data.table::data.table(matrix(rep(0, nids * nrow), nrow = nrow, ncol = nids))
  colnames(df) <- unlist(id_map)

  # All singletons.
  obs <- gDRutils:::.modify_polymapped_identifiers(df, exp_one_ids, id_map)
  expect_equal(obs, id_map)

  # Some polymappings.
  some_poly_map <- id_map
  some_poly_map[[1]] <- c(some_poly_map[[1]], "extra_item") 
  obs <- gDRutils:::.modify_polymapped_identifiers(df, exp_one_ids, id_map = some_poly_map)
  expect_equal(obs, id_map)

  # All polymappings.
  all_poly_map <- id_map
  all_poly_map <- lapply(seq_along(all_poly_map), function(x) c(all_poly_map[[x]], "extra_item"))
  names(all_poly_map) <- names(id_map)
  obs <- gDRutils:::.modify_polymapped_identifiers(df, exp_one_ids, id_map = all_poly_map)
  expect_equal(obs, id_map)

  # Multiple valid polymappings.
  multi_valid_map <- id_map
  multi_valid_map[[1]] <- c(multi_valid_map[[1]], multi_valid_map[[2]])
  expect_error(gDRutils:::.modify_polymapped_identifiers(df, exp_one_ids, id_map = multi_valid_map),
    regex = sprintf("multiple valid identifier values found for identifier: '%s'", names(multi_valid_map)[1]))
})

test_that("validate_identifiers works as expected", {
  nrow <- 5
  exp_one_ids <- get_expect_one_identifiers()
  req_ids <- get_required_identifiers()
  nids <- length(exp_one_ids)

  ids <- as.list(LETTERS[seq(nids)])
  names(ids) <- exp_one_ids

  df <- data.table::data.table(matrix(rep(0, nids * nrow), nrow = nrow, ncol = nids))
  colnames(df) <- unlist(ids)

  # Single mapping all valid.
  expect_equal(validate_identifiers(df, identifiers = ids, req_ids = req_ids, exp_one_ids = exp_one_ids), ids)

  # Single mapping but does not exist.
  single_map_ids <- ids
  single_map_ids[1] <- "Cinderella"
  expect_error(validate_identifiers(df, identifiers = single_map_ids, req_ids = req_ids, exp_one_ids = exp_one_ids))

  # Polymapping one exists.
  poly_map_ids <- ids
  poly_map_ids[[1]] <- c(poly_map_ids[[1]], "BOGUS")
  expect_equal(validate_identifiers(df, identifiers = poly_map_ids, req_ids = req_ids, exp_one_ids = exp_one_ids), ids)

  # Polymapping both exists.
  multi_valid_ids <- ids
  multi_valid_ids[[1]] <- c(multi_valid_ids[[1]], multi_valid_ids[[2]])
  expect_error(validate_identifiers(df, identifiers = multi_valid_ids, req_ids = req_ids, exp_one_ids = exp_one_ids),
    regex = sprintf("multiple valid identifier values found for identifier: '%s'", names(ids)[1]))
})
