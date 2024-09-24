context("standardize_MAE")

mapping_vector <- as.character(seq_len(5))
names(mapping_vector) <- names(iris)

test_that("rename_bumpy works as expected",  {
  bumpy <- BumpyMatrix::splitAsBumpyMatrix(iris, iris$Species, iris$Sepal.Length)
  bumpy_renamed <- rename_bumpy(bumpy, mapping_vector)
  expect_equal(BumpyMatrix::commonColnames(bumpy_renamed), unname(mapping_vector))
})

test_that("rename_DFrame works as expected",  {
  dframe <- S4Vectors::DataFrame(iris)
  dframe_renamed <- rename_DFrame(dframe, mapping_vector)
  expect_equal(names(dframe_renamed), unname(mapping_vector))
})

test_that("standardize_se works as expected",  {
  se_original <- get_synthetic_data("finalMAE_combo_matrix_small")[[1]]
  se <- se_original
  se@metadata$identifiers$drug <- "druuug"
  se@metadata$identifiers$concentration2 <- "dose 2"
  rowData(se) <- rename_DFrame(rowData(se), c("Gnumber" = "druuug"))
  assay(se, "RawTreated") <- rename_bumpy(assay(se, "RawTreated"), c("Concentration_2" = "dose 2"))
  se_standardized <- standardize_se(se)
  expect_equal(convert_se_assay_to_dt(se_standardized, "RawTreated"),
               convert_se_assay_to_dt(se_original, "RawTreated"))
})

test_that("standardize_se works as expected with default = FALSE",  {
  se_original <- get_synthetic_data("finalMAE_combo_matrix_small")[[1]]
  se <- se_original
  se@metadata$identifiers$drug <- "druuug"
  se@metadata$identifiers$concentration2 <- "dose 2"
  rowData(se) <- rename_DFrame(rowData(se), c("Gnumber" = "druuug"))
  assay(se, "RawTreated") <- rename_bumpy(assay(se, "RawTreated"), c("Concentration_2" = "dose 2"))
  se_standardized <- standardize_se(standardize_se(se), use_default = FALSE)
  expect_equal(convert_se_assay_to_dt(se_standardized, "RawTreated"),
               convert_se_assay_to_dt(se, "RawTreated"))
})


test_that("standardize_MAE works as expected",  {
  mae_original <- get_synthetic_data("finalMAE_combo_matrix_small")
  mae <- mae_original
  mae[[1]]@metadata$identifiers$drug <- "druuug"
  mae[[2]]@metadata$identifiers$drug <- "druuug"
  mae[[1]]@metadata$identifiers$concentration2 <- "dose 2"
  rowData(mae[[1]]) <- rename_DFrame(rowData(mae[[1]]), c("Gnumber" = "druuug"))
  rowData(mae[[2]]) <- rename_DFrame(rowData(mae[[2]]), c("Gnumber" = "druuug"))
  assay(mae[[1]], "RawTreated") <- rename_bumpy(assay(mae[[1]], "RawTreated"), c("Concentration_2" = "dose 2"))
  mae_standardized <- standardize_mae(mae)
  expect_equal(convert_mae_assay_to_dt(mae_standardized, "RawTreated"),
               convert_mae_assay_to_dt(mae_original, "RawTreated"))
})

test_that("standardize_MAE works with polymapped identifiers",  {
  mae_original <- get_synthetic_data("finalMAE_combo_matrix_small")
  mae <- mae_original
  mae[[1]]@metadata$identifiers$drug <- c("druuug", "Drug")
  mae[[2]]@metadata$identifiers$drug <- c("druuug", "Drug")
  mae[[1]]@metadata$identifiers$concentration2 <- "dose 2"
  rowData(mae[[1]]) <- rename_DFrame(rowData(mae[[1]]), c("Gnumber" = "druuug"))
  rowData(mae[[2]]) <- rename_DFrame(rowData(mae[[2]]), c("Gnumber" = "druuug"))
  assay(mae[[1]], "RawTreated") <- rename_bumpy(assay(mae[[1]], "RawTreated"), c("Concentration_2" = "dose 2"))
  mae_standardized <- standardize_mae(mae)
  expect_equal(convert_mae_assay_to_dt(mae_standardized, "RawTreated"),
               convert_mae_assay_to_dt(mae_original, "RawTreated"))
})

test_that("colData/rowData refinement functions work as expected",  {
  mae <- get_synthetic_data("finalMAE_combo_matrix_small")
  expect_true(inherits(refine_coldata(SummarizedExperiment::colData(mae[[1]]), mae[[1]]),
               "DataFrame"))
  expect_true(inherits(refine_rowdata(SummarizedExperiment::rowData(mae[[1]]), mae[[1]]),
               "DataFrame"))
  
  expect_error(refine_coldata(mae, mae), "Assertion on 'se' failed:")
  expect_error(refine_rowdata(mae, mae), "Assertion on 'se' failed:")
  
})


test_that("get_optional_rowdata_fields works as expected", {
  se <- get_synthetic_data("finalMAE_combo_matrix_small")[[1]]
  idfs <- get_SE_identifiers(se)
  opt_idfs <- get_optional_rowdata_fields(se)
  expect_equal(opt_idfs, unlist(idfs[c("drug_moa", "drug_moa2")],
                                use.names = FALSE))
  
  se2 <- get_synthetic_data("finalMAE_small")[[1]]
  idfs2 <- get_SE_identifiers(se2)
  opt_idfs2 <- get_optional_rowdata_fields(se2)
  expect_equal(opt_idfs2, idfs[["drug_moa"]])
})


test_that("set_unique_cl_names_dt and set_unique_drug_names_dt works correctly", {
  
  # DataFrame
  ## Duplicated CellLineName
  col_data <- S4Vectors::DataFrame(CellLineName = c("ID1", "ID1"), clid = c("C1", "C2"))
  res_1 <- set_unique_drug_names_dt(col_data)
  res_2 <- set_unique_cl_names_dt(col_data)
  expect_equal(col_data, res_1)
  expect_false(identical(col_data, res_2))
  expect_equal(c("ID1 (C1)", "ID1 (C2)"), res_2$CellLineName)
  
  ## Duplicated DrugName
  row_data <- S4Vectors::DataFrame(DrugName = c("DrugA", "DrugA"), Gnumber = c("G1", "G2"))
  res_3 <- set_unique_drug_names_dt(row_data)
  res_4 <- set_unique_cl_names_dt(row_data)
  expect_false(identical(row_data, res_3))
  expect_equal(row_data, res_4)
  expect_equal(c("DrugA (G1)", "DrugA (G2)"), res_3$DrugName)
  
  # data.table
  ## All different
  dt <- data.table::data.table(
    DrugName = c("DrugA", "DrugB", "DrugC", "DrugD"), 
    Gnumber = c("G1", "G2", "G3", "G4"),
    CellLineName = c("ID1", "ID2", "ID3", "ID4"), 
    clid = c("C1", "C2", "C3", "C4")
  )
  res_5 <- set_unique_drug_names_dt(dt)
  res_6 <- set_unique_cl_names_dt(dt)
  expect_equal(res_5, dt)
  expect_equal(res_6, dt)
  
  ## Duplicated CellLineName
  dt <- data.table::data.table(
    DrugName = c("DrugA", "DrugB", "DrugC", "DrugD", "DrugC", "DrugD"), 
    Gnumber = c("G1", "G2", "G3", "G4", "G3", "G4"),
    CellLineName = c("ID1", "ID1", "ID2", "ID2", "ID2", "ID2"), 
    clid = c("C1", "C2", "C3", "C4", "C5", "C6")
  )
  res_7 <- set_unique_drug_names_dt(dt)
  res_8 <- set_unique_cl_names_dt(dt)
  expect_equal(res_7, dt)
  expect_false(identical(res_8, dt))
  expect_equal(length(unique(res_8$CellLineName)), 5)

  ## Duplicated DrugName
  dt <- data.table::data.table(
    DrugName = c("DrugA", "DrugA", "DrugB", "DrugB", "DrugB", "DrugB"), 
    Gnumber = c("G1", "G2", "G3", "G4", "G5", "G6"),
    CellLineName = c("ID1", "ID2", "ID3", "ID4", "ID3", "ID4"), 
    clid = c("C1", "C2", "C3", "C4", "C3", "C4")
  )
  res_9 <- set_unique_drug_names_dt(dt)
  res_10 <- set_unique_cl_names_dt(dt)
  expect_false(identical(res_9, dt))
  expect_equal(length(unique(res_9$DrugName)), 6)
  expect_equal(res_10, dt)
  
  ## Duplicated both
  dt <- data.table::data.table(
    DrugName = c("DrugA", "DrugA", "DrugB", "DrugB"), 
    Gnumber = c("G1", "G2", "G3", "G4"),
    CellLineName = c("ID1", "ID1", "ID2", "ID2"), 
    clid = c("C1", "C2", "C3", "C4")
  )
  res_11 <- set_unique_drug_names_dt(dt)
  res_12 <- set_unique_cl_names_dt(dt)
  expect_false(identical(res_11, dt))
  expect_equal(length(unique(res_11$DrugName)), 4)
  expect_equal(length(unique(res_11$CellLineName)), 2)
  expect_false(identical(res_12, dt))
  expect_equal(length(unique(res_12$DrugName)), 2)
  expect_equal(length(unique(res_12$CellLineName)), 4)
  
})


test_that("set_unique_cl_names works correctly", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:4, ncol = 2)),
    colData = S4Vectors::DataFrame(CellLineName = c("ID1", "ID1"), clid = c("C1", "C2"))
  )
  se <- set_unique_cl_names(se)
  
  expect_equal(SummarizedExperiment::colData(se)$CellLineName, c("ID1 (C1)", "ID1 (C2)"))
})


test_that("set_unique_drug_names works correctly", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:4, ncol = 2)),
    rowData = S4Vectors::DataFrame(DrugName = c("DrugA", "DrugA"), Gnumber = c("G1", "G2"))
  )
  se <- set_unique_drug_names(se)
  
  expect_equal(SummarizedExperiment::rowData(se)$DrugName, c("DrugA (G1)", "DrugA (G2)"))
  
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:9, ncol = 3)),
    rowData = S4Vectors::DataFrame(DrugName = c("DrugA", "DrugA", "DrugB"),
                                   Gnumber = c("G1", "G2", "G5"),
                                   DrugName_2 = c("DrugC", "DrugC", "DrugD"),
                                   Gnumber_2 = c("G3", "G3", "G5")
    ))
  
  se2 <- set_unique_drug_names(se2)
  expect_equal(SummarizedExperiment::rowData(se2)$DrugName, c("DrugA (G1)", "DrugA (G2)", "DrugB"))
  expect_equal(SummarizedExperiment::rowData(se2)$DrugName_2, c("DrugC", "DrugC", "DrugD"))
})

test_that("set_unique_identifiers works correctly", {
  se1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:4, ncol = 2)),
    colData = S4Vectors::DataFrame(CellLineName = c("ID1", "ID1"), clid = c("C1", "C2")),
    rowData = S4Vectors::DataFrame(DrugName = c("DrugA", "DrugA"), Gnumber = c("G1", "G2"))
  )
  rownames(SummarizedExperiment::colData(se1)) <- c("Sample1", "Sample2")
  rownames(SummarizedExperiment::rowData(se1)) <- c("Gene1", "Gene2")
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(5:8, ncol = 2)),
    colData = S4Vectors::DataFrame(CellLineName = c("ID2", "ID2"), clid = c("C3", "C4")),
    rowData = S4Vectors::DataFrame(DrugName = c("DrugB", "DrugB"), Gnumber = c("G3", "G4"))
  )
  rownames(SummarizedExperiment::colData(se2)) <- c("Sample3", "Sample4")
  rownames(SummarizedExperiment::rowData(se2)) <- c("Gene3", "Gene4")
  mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(se1 = se1, se2 = se2))
  mae <- set_unique_identifiers(mae)
  
  expect_equal(SummarizedExperiment::colData(mae[[1]])$CellLineName, c("ID1 (C1)", "ID1 (C2)"))
  expect_equal(SummarizedExperiment::rowData(mae[[1]])$DrugName, c("DrugA (G1)", "DrugA (G2)"))
  expect_equal(SummarizedExperiment::colData(mae[[2]])$CellLineName, c("ID2 (C3)", "ID2 (C4)"))
  expect_equal(SummarizedExperiment::rowData(mae[[2]])$DrugName, c("DrugB (G3)", "DrugB (G4)"))
})
