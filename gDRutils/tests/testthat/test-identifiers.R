#library(testthat); library(gDRutils)

context("test functions relating to getting, setting, and resetting headers")

test_that("get_identifier and set_identifier work", {
  reset_identifiers()

  expect_error(get_identifier("BOGUS"))

  expect_equal(get_identifier("duration"), "Duration")

  set_identifier("cellline", "my_personal_cell_line_identifiers")
  expect_equal(get_identifier("cellline"), "my_personal_cell_line_identifiers")

  expect_equal(get_identifier(), list(duration = "Duration", 
                                      cellline = "my_personal_cell_line_identifiers", 
                                      cellline_name = "CellLineName",
                                      cellline_tissue = "Tissue",
                                      cellline_ref_div_time = "ReferenceDivisionTime",
                                      drug = "Gnumber",
                                      drugname = "DrugName",
                                      untreated_tag = c("untreated", "vehicle"), 
                                      masked_tag = "masked", 
                                      WellPosition = c("WellRow", "WellColumn"), 
                                      well_position = c("WellRow", "WellColumn")))
})


test_that("reset_identifiers works", {
  reset_identifiers()

  d <- get_identifier()
  expect_true(length(d) != 0L)
  set_identifier("duration", "TEST_DURATION")
  expect_equal(get_identifier("duration"), "TEST_DURATION")

  ilist <- reset_identifiers()
  expect_equal(ilist, NULL)
  expect_equal(get_identifier("duration"), "Duration")
})
