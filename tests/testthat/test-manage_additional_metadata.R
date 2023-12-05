context("Tests for helpers-manageData")
m <- 9
n <- 5

synthetic_data <- gen_synthetic_data(m, n)
drug_names <- synthetic_data$drug_names
cell_names <- synthetic_data$cell_names
dt <- synthetic_data$dt

test_that("addClass works as expected", {
  expect_s3_class(addClass(dt, "testingClass"), "testingClass")
})
