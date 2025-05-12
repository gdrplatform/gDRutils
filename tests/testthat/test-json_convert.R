library(testthat)

test_that(".standardize_column_names works as expected", {
  ids <- letters
  names(ids) <- LETTERS
  df <- data.table::data.table(b = seq(100),
                             m = seq(100),
                             a = seq(100),
                             aa = seq(100))
  obs <- .standardize_column_names(df, ids[c("A", "B", "M")])
  expect_equal(colnames(obs), c("B", "M", "A", "aa"))
})

test_that(".convert_element_metadata_to_json works as expected", {
  m <- 100
  n <- 26
  data <- data.table::as.data.table(matrix(rep(seq(m), n), nrow = m))
  shuffle <- sample(1:n, n)
  data.table::setnames(data, LETTERS[shuffle])

  expect_error(.convert_element_metadata_to_json(data, req_cols = c("AA", "A", "B")))

  req_exp <- LETTERS[1:5]
  obs <- .convert_element_metadata_to_json(data, req_cols = req_exp)
  obs_main_json <- jsonlite::fromJSON(sprintf("{%s}", obs$main))
  expect_equal(names(obs_main_json), req_exp)

  obs_opt_json <- jsonlite::fromJSON(obs$opt)
  obs_opt <- setdiff(LETTERS, names(obs_opt_json))
  expect_equal(obs_opt, req_exp)
})

test_that("convert_se_to_json works with diferent types of metadata", {
  rdata <-
    data.table::data.table(
      mydrug = letters,
      mydrugname = letters,
      mydrugmoa = letters,
      Duration = 1
    )
  cdata <-
    data.table::data.table(
      mycellline = letters,
      mycelllinename = letters,
      mycelllinetissue = letters,
      cellline_ref_div_time = letters
    )
  identifiers <- list(cellline = "mycellline",
                      cellline_name = "mycelllinename",
                      cellline_tissue = "mycelllinetissue",
                      cellline_ref_div_time = "cellline_ref_div_time",
                      drug = "mydrug",
                      drug_name = "mydrugname",
                      drug_moa = "mydrugmoa",
                      duration = "Duration")
  se <- SummarizedExperiment::SummarizedExperiment(rowData = rdata,
                                                   colData = cdata)
  se <- set_SE_identifiers(se, identifiers)

  # list
  md <- list(title = "my awesome experiment",
             description = "description of experiment",
             source = list(name = "GeneData_Screener", id = "QCS-12345"))
  se <- set_SE_experiment_metadata(se, md)
  # check if there is no error
  expect_error(convert_se_to_json(se), NA)

  # data.table
  md <- data.table::data.table(title = "my awesome experiment",
                       description = "description of experiment",
                       source = "GeneData_Screener")
  se <- set_SE_experiment_metadata(se, md)

  # check if there is no error
  expect_error(convert_se_to_json(se), NA)
})

test_that("strip_first_and_last_char works as expected", {
  expect_equal(strip_first_and_last_char("hello"), "ell")
  expect_equal(strip_first_and_last_char("{}"), "")
})

