
context("Tests for helpers-manageData")
m <- 9
n <- 5

synthetic_data <- gen_synthetic_data(m, n)
drug_names <- synthetic_data$drug_names
cell_names <- synthetic_data$cell_names
dt <- synthetic_data$dt

test_that("addClass works as expected", {
  expect_class(addClass(dt, "testingClass"), "testingClass")
})


test_that("modifyData works as expected", {
  dt2 <- dt
  dt2[["GR_AOC (GDS)"]] <- 1
  data <- list(dt2)
  dtReplaced <- replaceMetrics(list(dt2), c("gDR", "GDS"), "GDS")
  expect_subset("GR_AOC (gDR)", names(dtReplaced[[1]]))
  expect_equal(dt2[["GR_AOC"]], dtReplaced[[1]][["GR_AOC (gDR)"]])
  expect_equal(dt2[["GR_AOC (GDS)"]], dtReplaced[[1]][["GR_AOC"]])
  
  dtRevert <- replaceMetrics(dtReplaced, c("gDR", "GDS"), "gDR")
  expect_identical(dtRevert[[1]][, sort(names(dtRevert[[1]]))], dt2[, sort(names(dt2))])
})
