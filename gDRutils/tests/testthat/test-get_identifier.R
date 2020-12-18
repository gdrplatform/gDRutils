#library(testthat); library(gDRutils)

test_that("get_identifier works", {
  expect_error(get_identifier("BOGUS"))

  expect_equal(get_identifier("duration"), "Duration")
  expect_error(set_identifier("duration", "BLAH"))

  set_identifier("cellline", "my_personal_cell_line_identifiers")
  expect_equal(get_identifier("cellline"), "my_personal_cell_line_identifiers")

  expect_equal(get_identifier(), list(duration = "Duration", 
                                      cellline = "my_personal_cell_line_identifiers", 
                                      drug = "Gnumber",
                                      drugname = "DrugName",
                                      untreated_tag = c("untreated", "vehicle"), 
                                      masked_tag = "masked", 
                                      WellPosition = c("WellRow", "WellColumn")))
})
