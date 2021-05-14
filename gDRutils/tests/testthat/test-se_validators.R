test_that("validate_se_assay_name works as expected", {
  se <- SummarizedExperiment::SummarizedExperiment(assays = list("orange" = matrix(1, 1, 1)))
  expect_error(validate_se_assay_name(se, "apple"))
  expect_true(is.null(validate_se_assay_name(se, "orange")))
})

test_that("validate_se works as expected", {
  se1 <- SummarizedExperiment::SummarizedExperiment(assays = list("orange" = matrix(1, 1, 1)))
  expect_error(validate_SE(se1))
  x <- IRanges::NumericList(split(runif(1000), factor(sample(50, 1000, replace = TRUE), 1:50))) 
  se2 <- SummarizedExperiment::SummarizedExperiment(assays = 
                                                      list("RawTreated" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "Controls" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "Normalized" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "RefGRvalue" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "RefRelativeViability" = 
                                                             BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "DivisionTime" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "Averaged" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "Metrics" = BumpyMatrix::BumpyMatrix(x, c(10, 5))))
  
  expect_error(validate_SE(se2))
  S4Vectors::metadata(se2) <- vector(mode = "list", length = 7)
  names(S4Vectors::metadata(se2)) <- c("experiment_metadata",
                                       "df_",
                                       "Keys", "df_raw_data",
                                       "fit_parameters",
                                       "drug_combinations",
                                       ".internals")
  expect_error(validate_SE(se2))
  
})
