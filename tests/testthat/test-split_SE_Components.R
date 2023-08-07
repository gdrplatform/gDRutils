test_that("split_SE_components splits the correct columns", {
  # Standard case.
  md <- split_SE_components(test_df)
  expect_true(all(c("Gnumber", "DrugName", "replicates", "drug_moa") %in% colnames(md$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% colnames(md$condition_md)))
  expect_equal(sum(ncol(md$treatment_md), ncol(md$condition_md), length(md$data_fields), ncol(md$experiment_md)), 
    ncol(test_df))
  pure <- get_env_identifiers(simplify = TRUE)
  expect_equal(md$identifiers_md[names(pure)], pure)

  # nested_keys argument works as expected.
  md2 <- split_SE_components(test_df, nested_keys = c("replicates"), combine_on = 1)
  expect_true(all(c("Gnumber", "DrugName", "drug_moa") %in% colnames(md2$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% colnames(md2$condition_md)))
  expect_true(all(c("WellRow", "WellColumn", "replicates") %in% md2$data_fields))
  expect_equal(ncol(test_df), 
    sum(ncol(md2$treatment_md), ncol(md2$condition_md), length(md2$data_fields), ncol(md2$experiment_md)))

  # combine_on argument works as expected
  md3 <- split_SE_components(test_df, nested_keys = c("replicates"), combine_on = 2L)
  expect_true(all(c("Gnumber", "DrugName", "drug_moa") %in% colnames(md3$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% 
    colnames(md3$condition_md)))
  expect_true(all(c("WellRow", "WellColumn", "replicates") %in% md3$data_fields))
  expect_equal(ncol(test_df), 
    sum(ncol(md3$treatment_md), ncol(md3$condition_md), length(md3$data_fields), ncol(md3$experiment_md)))

  # nested key is a main identifier.
  md4 <- split_SE_components(test_df, nested_keys = c("drug_moa"))
  expect_true(all(c("Gnumber", "DrugName", "replicates") %in% colnames(md4$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% colnames(md4$condition_md)))
  expect_true(all(c("WellRow", "WellColumn", "drug_moa") %in% md4$data_fields))
  expect_equal(ncol(test_df), 
    sum(ncol(md4$treatment_md), ncol(md4$condition_md), length(md4$data_fields), ncol(md4$experiment_md)))
  
  # order of columns is correct
  expect_equal(names(md$treatment_md), c("Gnumber", "DrugName", "drug_moa", "Duration", "replicates"))
  
  # split_SE_components with changed identifiers
  new_identifier_name <- "SomeDrug"
  set_env_identifier("drug", new_identifier_name)
  test_df_modified <- test_df
  names(test_df_modified)[1] <- new_identifier_name
  md5 <- split_SE_components(test_df_modified)
  expect_equal(names(md5$treatment_md), c(new_identifier_name, "DrugName", "drug_moa", "Duration", "replicates"))
})


test_that("split_SE_components throws a warning for bad cell line metadata", {
  df2 <- test_df
  df2$CellLineName[7] <- "Abnormality"
  expect_warning(split_SE_components(df2), 
    regexp = "'CellLineName' not metadata for unique cell line identifier column")
})


test_that("add_rownames_to_metadata works as expected", {
  cols <- c("a", "b")
  md <- data.frame(a = LETTERS, b = letters, c = paste0(LETTERS, letters))
  expect_true(all(rownames(md) == as.character(seq(nrow(md)))))
  out <- gDRutils:::add_rownames_to_metadata(md, cols)
  expect_true(all(rownames(out) != as.character(seq(nrow(md)))))
  expect_equal(colnames(out), cols)
})

test_that("split_SE_components returns rowData in a proper order", {
  md <- split_SE_components(test_df)
  expect_true(all(gsub("_.*", "",
                       rownames(md$treatment_md)) ==
                    md$treatment_md[[gDRutils::get_env_identifiers("drug")]]))
})

test_that("split_SE_components works with colnames with -", {
  test_df2 <- data.table::copy(test_df)
  test_df2$`fix5-aza` <- sample(c(0.5 ,0), size = nrow(test_df2), replace = TRUE)
  md <- split_SE_components(test_df2)
  expect_true("fix5-aza" %in% names(md$treatment_md))
})


