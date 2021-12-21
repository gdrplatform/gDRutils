library(testthat)
context("assay_name")

test_that("get_assay_names",  {
  ### without SE
  ## all values
  agan <- get_assay_names()
  expect_true(length(agan) > 2)
  expect_identical("character", class(agan))
  expect_true(!is.null(names(agan)))
  
  ## subset
  # single value 
  sgan <- get_assay_names(type = names(agan[1]))
  expect_identical(agan[[1]],sgan)
  sgan <- get_assay_names(type = names(agan[1]), simplify = FALSE)
  expect_identical(agan[1],sgan)
  sgan <- get_assay_names(type = names(agan[1]), prettify = TRUE)
  expect_false(agan[1] == sgan)
  # multiple values
  sgan <- get_assay_names(type = names(agan[1:2]))
  expect_identical(agan[1:2],sgan)
  
 
  ### errors 
  # bad value for given filter provided
  expect_error(get_assay_names(type = "bad_type"), "Assertion on 'bad_type' failed")
  expect_error(get_assay_names(group = "inv_group"), "Assertion on 'inv_group' failed")
  expect_error(get_assay_names(data_type = "inv_data_type"), "Assertion on 'inv_data_type' failed")

})
