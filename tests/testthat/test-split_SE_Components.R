test_that("split_SE_components splits the correct columns", {
  # Standard case.
  md <- split_SE_components(test_df)
  expect_true(all(c("Gnumber", "DrugName", "drug_moa") %in% colnames(md$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "Replicate", "ReferenceDivisionTime") %in%
                    colnames(md$condition_md)))
  expect_equal(sum(NCOL(md$treatment_md), NCOL(md$condition_md), length(md$data_fields), NCOL(md$experiment_md)),
    NCOL(test_df))
  pure <- get_env_identifiers(simplify = TRUE)
  expect_equal(md$identifiers_md[names(pure)], pure)

  # nested_keys argument works as expected.
  md2 <- split_SE_components(test_df, nested_keys = "Replicate", combine_on = 1)
  expect_true(all(c("Gnumber", "DrugName", "drug_moa") %in% colnames(md2$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% colnames(md2$condition_md)))
  expect_true(all(c("WellRow", "WellColumn", "Replicate") %in% md2$data_fields))
  expect_equal(NCOL(test_df),
    sum(NCOL(md2$treatment_md), NCOL(md2$condition_md), length(md2$data_fields), NCOL(md2$experiment_md)))

  # combine_on argument works as expected
  md3 <- split_SE_components(test_df, nested_keys = "Replicate", combine_on = 2L)
  expect_true(all(c("Gnumber", "DrugName", "drug_moa") %in% colnames(md3$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in%
    colnames(md3$condition_md)))
  expect_true(all(c("WellRow", "WellColumn", "Replicate") %in% md3$data_fields))
  expect_equal(NCOL(test_df),
    sum(NCOL(md3$treatment_md), NCOL(md3$condition_md), length(md3$data_fields), NCOL(md3$experiment_md)))

  # nested key is a main identifier.
  md4 <- split_SE_components(test_df, nested_keys = c("drug_moa"))
  expect_true(all(c("Gnumber", "DrugName") %in% colnames(md4$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Replicate", "Tissue", "ReferenceDivisionTime") %in%
                    colnames(md4$condition_md)))
  expect_true(all(c("WellRow", "WellColumn", "drug_moa") %in% md4$data_fields))
  expect_equal(NCOL(test_df),
    sum(NCOL(md4$treatment_md), NCOL(md4$condition_md), length(md4$data_fields), NCOL(md4$experiment_md)))

  # order of columns is correct
  expect_equal(names(md$treatment_md), c("Gnumber", "DrugName", "drug_moa", "Duration"))

  # split_SE_components with changed identifiers
  new_identifier_name <- "SomeDrug"
  set_env_identifier("drug", new_identifier_name)
  test_df_modified <- test_df
  names(test_df_modified)[1] <- new_identifier_name
  md5 <- split_SE_components(test_df_modified)
  expect_equal(names(md5$treatment_md), c(new_identifier_name, "DrugName", "drug_moa", "Duration"))
})

test_that("add_rownames_to_metadata works as expected", {
  cols <- c("a", "b")
  md <- data.frame(a = LETTERS, b = letters, c = paste0(LETTERS, letters))
  expect_true(all(rownames(md) == as.character(seq_len(NROW(md)))))
  out <- add_rownames_to_metadata(md, cols)
  expect_true(all(rownames(out) != as.character(seq_len(NROW(md)))))
  expect_equal(colnames(out), cols)
})

test_that("split_SE_components returns rowData in a proper order", {
  md <- split_SE_components(test_df)
  expect_true(all(gsub("_.*", "",
                       rownames(md$treatment_md)) ==
                    md$treatment_md[[get_env_identifiers("drug")]]))
})

test_that("split_SE_components works with colnames with -", {
  test_df2 <- data.table::copy(test_df)
  test_df2$`fix5-aza` <- sample(c(0.5, 0), size = NROW(test_df2), replace = TRUE)
  md <- split_SE_components(test_df2)
  expect_true("fix5-aza" %in% names(md$treatment_md))
})

test_that("split_SE_components always returns experiment_md as a list", {
  # Case 1: no constant columns -> empty list
  df_no_const <- data.table::data.table(
    clid = c("CL1", "CL2"),
    Gnumber = c("DrugA", "DrugB"),
    ReadoutValue = c(1.0, 2.0)
  )
  md_no_const <- split_SE_components(df_no_const)
  expect_true(is.list(md_no_const$experiment_md))
  expect_false(is.data.frame(md_no_const$experiment_md))
  expect_equal(length(md_no_const$experiment_md), 0L)

  # Case 2: constant columns present -> named list with the constant value
  # BatchNumber is not a known identifier so it lands in experiment_md when constant
  df_const <- data.table::data.table(
    clid = c("CL1", "CL2"),
    Gnumber = c("DrugA", "DrugB"),
    Duration = c(72, 72),
    BatchNumber = c("B1", "B1"),
    ReadoutValue = c(1.0, 2.0)
  )
  md_const <- split_SE_components(df_const)
  expect_true(is.list(md_const$experiment_md))
  expect_false(is.data.frame(md_const$experiment_md))
  expect_true(length(md_const$experiment_md) > 0L)
  expect_equal(md_const$experiment_md[["BatchNumber"]], "B1")

  # Case 3: standard test_df -> still a list
  md_std <- split_SE_components(test_df)
  expect_true(is.list(md_std$experiment_md))
  expect_false(is.data.frame(md_std$experiment_md))
})

test_that("split_SE_components sorts non-default columns", {
  test_df3 <- data.table::copy(test_df)
  test_df3$`fix5-aza` <- sample(c(0.5, 0), size = NROW(test_df3), replace = TRUE)
  test_df3$`a-fix5-aza` <- sample(c(0.5, 0), size = NROW(test_df3), replace = TRUE)
  test_df3$`b-fix5-aza` <- sample(c(0.5, 0), size = NROW(test_df3), replace = TRUE)
  md <- split_SE_components(test_df3)
  expect_identical(grep("fix5-aza", names(md$treatment_md), value = TRUE),
                   c("a-fix5-aza", "b-fix5-aza", "fix5-aza"))

  md2 <- split_SE_components(test_df3, combine_on = 2)
  expect_identical(grep("fix5-aza", names(md2$condition_md), value = TRUE),
                   c("a-fix5-aza", "b-fix5-aza", "fix5-aza"))

  # ensure identical rownames of DFrames
  # in case of two data.tables with assay data with different order of non-default columns
  test_df4 <- data.table::copy(test_df3)
  data.table::setcolorder(test_df4, c("b-fix5-aza", "a-fix5-aza", "fix5-aza"))
  md3 <- split_SE_components(test_df4)
  expect_identical(sort(rownames(md$treatment_md)), sort(rownames(md3$treatment_md)))
})
