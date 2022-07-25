library(testthat)

test_that("validate works as expected", {
  identifiers <- list(drug = "drug", drug_name = "drug_name", drug_moa = "drug_moa", 
                       cellline = "cellline", cellline_name = "cellline_name", cellline_tissue = "cellline_tissue")
  rdata <- S4Vectors::DataFrame(drug = c("A", "B"), 
                                drug_name = c("AA", "AB"),
                                drug_moa = c("action1", "action2"),
                                extra_row = c("Bruce", "Lee"),
                                row.names = seq(2))

  cdata <- S4Vectors::DataFrame(cellline = c("ID1", "ID2"),
                                cellline_name = c("ABC", "DEF"),
                                cellline_tissue = c("Breast", "Lung"),
                                extra_col = c("Barack", "Obama"),
                                row.names = seq(2))

  se <- SummarizedExperiment::SummarizedExperiment(rowData = rdata, colData = cdata)
  se <- gDRutils::set_SE_identifiers(se, identifiers)
  md <- list(sources = list(list(name = "Screen Gene Data", id = "ID-12345"),
             list(name = "University_of_Genes", id = "UT-1234")),
             description = "description test",
             title = "title test",
             experimentalist = "abcde",
             duration = 4)
  se <- gDRutils::set_SE_experiment_metadata(se, md)
  sejson <- convert_se_to_json(se)
  schema_path <- file.path(system.file("schemas", package = "gDRelastic"), "se.json")
  expect_true(validate_json(sejson, schema_path = schema_path))
  
  # original MAE not converted to DSDB MAE
  # errors expected
  tmae1 <-
    readRDS(system.file(package = "gDRtestData", "testdata/finalMAE_small.RDS"))
  v_st <- validate_mae_with_schema(tmae1)
  exp_names_v <- paste0("experiment:", names(MultiAssayExperiment::experiments(tmae1)))
  expect_true(all(exp_names_v %in% names(v_st)))
  expect_false(v_st[[1]])
  
  v_st_att <- attributes(v_st[[1]])
  expect_identical(v_st_att$exit_code, 2)
  expect_true(nchar(v_st_att$error) > 0)
  expect_true(inherits(v_st_att$derror, "data.frame"))
  expect_true(nrow(v_st_att$derror) > 0)
 
  # fix validation errors by adding proper metadata 
  S4Vectors::metadata(tmae1[[1]]) <-
    list(
      identifiers = get_SE_identifiers(tmae1[[1]]),
      experiment_metadata = list(
        title = "tt",
        description = "dummy_desc",
        experimentalist = "doej",
        duration = 10,
        sources = list(list(id = "testId", name = "dummyName"))
      )
    )
  v_st <- validate_mae_with_schema(tmae1)
  expect_true(v_st[[1]])
  
})
