library(testthat)

test_that(".standardize_column_names works as expected", {
  ids <- letters
  names(ids) <- LETTERS
  df <- S4Vectors::DataFrame(b = seq(100), 
                             m = seq(100),
                             a = seq(100),
                             aa = seq(100))
  obs <- .standardize_column_names(df, ids[c("A", "B", "M")])
  expect_equal(colnames(obs), c("B", "M", "A", "aa"))
})

test_that(".convert_element_metadata_to_json works as expected", {
  m <- 100
  n <- 26
  data <-  as.data.frame(matrix(rep(seq(m), n), nrow = m))
  shuffle <- sample(1:n, n)
  names(data) <- LETTERS[shuffle]

  expect_error(.convert_element_metadata_to_json(data, req_cols = c("AA", "A", "B")))

  req_exp <- LETTERS[1:5]
  obs <- .convert_element_metadata_to_json(data, req_cols = req_exp)
  obs_main_json <- jsonlite::fromJSON(sprintf("{%s}", obs$main))
  expect_equal(names(obs_main_json), req_exp)

  obs_opt_json <- jsonlite::fromJSON(obs$opt)
  obs_opt <- setdiff(LETTERS, names(obs_opt_json))
  expect_equal(obs_opt, req_exp)
})

test_that("strip_first_and_last_char works as expected", {
  expect_equal(strip_first_and_last_char("hello"), "ell")
})
