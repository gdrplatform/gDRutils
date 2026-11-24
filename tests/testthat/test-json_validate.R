library(testthat)

test_that("validate works as expected", {
  identifiers <- list(drug = "drug", drug_name = "drug_name", drug_moa = "drug_moa", duration = "duration",
                       cellline = "cellline", cellline_name = "cellline_name", cellline_tissue = "cellline_tissue",
                      cellline_ref_div_time = "cellline_ref_div_time")
  rdata <- data.table::data.table(drug = c("A", "B"),
                                drug_name = c("AA", "AB"),
                                drug_moa = c("action1", "action2"),
                                duration = 1:2,
                                extra_row = c("Bruce", "Lee"),
                                row.names = seq(2))

  cdata <- data.table::data.table(cellline = c("ID1", "ID2"),
                                cellline_name = c("ABC", "DEF"),
                                cellline_tissue = c("Breast", "Lung"),
                                cellline_ref_div_time = 1:2,
                                extra_col = c("Barack", "Obama"),
                                row.names = seq(2))

  se <- SummarizedExperiment::SummarizedExperiment(rowData = rdata, colData = cdata)
  se <- set_SE_identifiers(se, identifiers)
  md <- list(sources = list(list(name = "Screen Gene Data", id = "ID-12345"),
             list(name = "University_of_Genes", id = "UT-1234")),
             description = "description test",
             title = "title test",
             experimentalist = "abcde")
  se <- set_SE_experiment_metadata(se, md)
  sejson <- convert_se_to_json(se)
  schema_path <- file.path(system.file("schemas", package = "gDRutils"), "se.json")
  expect_true(validate_json(sejson, schema_path = schema_path))

  # original MAE not converted to DSDB MAE
  # errors expected
  tmae1 <-
    get_synthetic_data("finalMAE_small")
  v_st <- validate_mae_with_schema(tmae1)
  exp_names_v <- c("mae", paste0("experiment:", names(MultiAssayExperiment::experiments(tmae1))))
  expect_true(all(exp_names_v %in% names(v_st)))
  expect_true(v_st[["mae"]])
  expect_false(v_st[[exp_names_v[2]]])

  v_st_att <- attributes(v_st[[exp_names_v[2]]])
  expect_identical(v_st_att$exit_code, 2)
  expect_true(nchar(v_st_att$error) > 0)
  expect_true(inherits(v_st_att$derror, "data.frame"))
  expect_true(NROW(v_st_att$derror) > 0)

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
  expect_true(v_st[[exp_names_v[2]]])

})

test_that("validate accepts an unknown duration", {
  identifiers <- list(drug = "drug", drug_name = "drug_name", drug_moa = "drug_moa", duration = "duration",
                      cellline = "cellline", cellline_name = "cellline_name", cellline_tissue = "cellline_tissue",
                      cellline_ref_div_time = "cellline_ref_div_time")
  make_se <- function(duration) {
    rdata <- data.table::data.table(drug = c("A", "B"), drug_name = c("AA", "AB"),
                                    drug_moa = c("action1", "action2"), duration = duration)
    cdata <- data.table::data.table(cellline = c("ID1", "ID2"), cellline_name = c("ABC", "DEF"),
                                    cellline_tissue = c("Breast", "Lung"), cellline_ref_div_time = c(NA, 2))
    se <- SummarizedExperiment::SummarizedExperiment(rowData = rdata, colData = cdata)
    se <- set_SE_identifiers(se, identifiers)
    set_SE_experiment_metadata(se, list(sources = list(list(name = "PRISM", id = "ID-1")),
                                        description = "d", title = "t", experimentalist = "abcde"))
  }
  schema_path <- file.path(system.file("schemas", package = "gDRutils"), "se.json")

  # PRISM data has no assay time; NA reaches the JSON as the string "NA"
  expect_true(validate_json(convert_se_to_json(make_se(c(NA_real_, NA_real_))), schema_path = schema_path))
  expect_true(validate_json(convert_se_to_json(make_se(c(72, NA_real_))), schema_path = schema_path))
  # any other text is still rejected
  expect_false(validate_json(convert_se_to_json(make_se(c("72h", "NA"))), schema_path = schema_path))
})
