test_that("validate_se_assay_name works as expected", {
  se <- SummarizedExperiment::SummarizedExperiment(assays = list("orange" = matrix(1, 1, 1)))
  expect_error(validate_se_assay_name(se, "apple"))
  expect_true(is.null(validate_se_assay_name(se, "orange")))
})
