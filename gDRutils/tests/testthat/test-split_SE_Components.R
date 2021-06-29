library(testthat); library(gDRutils); source("setup.R")
test_that("split_SE_components splits the correct columns", {
  # Standard case.
  md <- split_SE_components(test_df)
  expect_true(all(c("Gnumber", "DrugName", "replicates", "drug_moa") %in% colnames(md$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% colnames(md$condition_md)))
  expect_equal(sum(ncol(md$treatment_md), ncol(md$condition_md), length(md$data_fields), ncol(md$experiment_md)), 
    ncol(test_df))
  pure <- get_identifier()
  expect_equal(md$identifiers_md[names(pure)], pure)

  # Check that nested_keys argument works as expected.
  md2 <- split_SE_components(test_df, nested_keys = c("replicates"))
  expect_true(all(c("Gnumber", "DrugName", "drug_moa", "Concentration") %in% colnames(md2$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% colnames(md2$condition_md)))
  expect_true(all(c("WellRow", "WellColumn", "replicates") %in% md2$data_fields))
  expect_equal(ncol(test_df), 
    sum(ncol(md2$treatment_md), ncol(md2$condition_md), length(md2$data_fields), ncol(md2$experiment_md)))
})


test_that("split_SE_components throws a warning for bad cell line metadata", {
  df2 <- test_df
  df2$CellLineName[7] <- "Abnormality"
  expect_warning(split_SE_components(df2), 
    regexp = "'CellLineName' not metadata for unique cell line identifier column")
})
