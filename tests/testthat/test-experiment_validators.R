test_that("validate_dimnames works as expected", {
  err_msg1 <-
    "Assertion on 'identical(dn1[idx], dn2[idx])' failed: Must be TRUE."
  expect_error(validate_dimnames(data.frame(a = 2, b = 3), data.frame(a = 3, c = 4)), err_msg1, fixed = TRUE)
  expect_error(validate_dimnames(data.frame(a = 2, b = 3), data.frame(a = 3, b = 4, c = 4)), err_msg1, fixed = TRUE)
  expect_error(validate_dimnames(data.frame(), data.frame(a = 3, b = 4, c = 4)), err_msg1, fixed = TRUE)
  
  expect_null(validate_dimnames(data.frame(a = 2, b = 3), data.frame(a = 3, b = 4)))
  expect_null(validate_dimnames(data.frame(), data.frame()))
})

test_that("validate_se_assay_name works as expected", {
  se <- SummarizedExperiment::SummarizedExperiment(assays = list("orange" = matrix(1, 1, 1)))
  expect_error(validate_se_assay_name(se, "apple"))
  expect_true(is.null(validate_se_assay_name(se, "orange")))
})

test_that("validate_SE works as expected", {
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


test_that("validate_SE works as expected on real data", {
  maeReal <- get_synthetic_data("finalMAE_small")

  validate_SE(maeReal[[1]])
  # Add empty drug_moa in one record of rowData
  rowdata <- SummarizedExperiment::rowData(maeReal[[1]])
  rowdata[2, "drug_moa"] <- ""
  SummarizedExperiment::rowData(maeReal[[1]]) <- rowdata
  expect_error(validate_SE(maeReal[[1]]),
               regexp = "Assertion on \'any(stats::na.omit(unlist(rowdata)) == \"\")\' failed: Must be FALSE.",
               fixed = TRUE)
})

test_that("validate_mae works as expected", {
  m <- 20
  n <- 10
  rnames <- LETTERS[1:m]
  cnames <- letters[1:n]
  
  ref_gr_value <-  matrix(runif(m * n), nrow = m, ncol = n, dimnames = list(rnames, cnames))
  se1 <- SummarizedExperiment::SummarizedExperiment(assays = list(RefGRvalue = ref_gr_value[1:10, ]),
                                                    rowData = S4Vectors::DataFrame(rnames)[1:10, , drop = FALSE],
                                                    colData = S4Vectors::DataFrame(cnames))
  
  se2 <- SummarizedExperiment::SummarizedExperiment(assays = list(RefGRvalue = ref_gr_value[11:20, ]),
                                                    rowData = S4Vectors::DataFrame(rnames)[11:20, , drop = FALSE],
                                                    colData = S4Vectors::DataFrame(cnames))
  
  mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(one = se1,
                                                                       two = se2))
  
  expect_error(validate_MAE(mae))
  x <- IRanges::NumericList(split(runif(1000), factor(sample(50, 1000, replace = TRUE), 1:50))) 
  se3 <- SummarizedExperiment::SummarizedExperiment(assays = 
                                                      list("RawTreated" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "Controls" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "Normalized" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "RefGRvalue" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "RefRelativeViability" = 
                                                             BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "DivisionTime" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "Averaged" = BumpyMatrix::BumpyMatrix(x, c(10, 5)),
                                                           "Metrics" = BumpyMatrix::BumpyMatrix(x, c(10, 5))))
  
  colData(se3) <- methods::new(
    "DFrame",
    rownames = c("A", "B", "C", "D", "E"),
    nrows = 5L,
    listData = list(Treatment = c("ChIP", "Input", "ChIP", "Input", "ChIP")),
    elementType = "ANY",
    elementMetadata = NULL,
    metadata = list()
  )
  
  mae2 <- MultiAssayExperiment::MultiAssayExperiment(experiments = list("single-agent" = se3))
  
  maeReal <- get_synthetic_data("finalMAE_small")
  validate_MAE(maeReal)
  maeReal2 <- MultiAssayExperiment::MultiAssayExperiment(experiments = list("single-agent" = maeReal[[1]],
                                                                            "matrix" = maeReal[[1]]))
  validate_MAE(maeReal2)
})
