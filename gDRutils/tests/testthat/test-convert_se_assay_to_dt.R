library(testthat)
library(gDRutils)

test_that("convert_se_assay_to_dt works as expected", {
  m <- 20
  n <- 10
  rnames <- LETTERS[1:m]
  cnames <- letters[1:n]
  
  # Normal matrix.
  ref_gr_value <-  matrix(runif(m * n), nrow = m, ncol = n, dimnames = list(rnames, cnames))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(RefGRvalue = ref_gr_value),
                                                   rowData = S4Vectors::DataFrame(rnames),
                                                   colData = S4Vectors::DataFrame(cnames))
  
  dt <- convert_se_assay_to_dt(se = se, assay_name = "RefGRvalue", include_metadata = FALSE)
  expect_equal(dt$RefGRvalue, as.vector(ref_gr_value))
  expect_equal(dim(dt), c(m * n, 3))
  
  dt <- convert_se_assay_to_dt(se = se, assay_name = "RefGRvalue", include_metadata = TRUE)
  expect_equal(dt[order(dt$cId), "RefGRvalue"][[1]], as.vector(ref_gr_value))
  expect_equal(dim(dt), c(m * n, 5))
  expect_equal(dt$rnames, as.character(dt$rId))
  expect_equal(dt$cnames, as.character(dt$cId))
  
  # BumpyDataFrameMatrix.
  df <- S4Vectors::DataFrame(r = rep(rnames, n), c = rep(cnames, m), values = runif(m * n))
  nested_rnames <- seq(nrow(df))
  rownames(df) <- nested_rnames
  
  norm <- BumpyMatrix::splitAsBumpyMatrix(df, row = df$r, column = df$c)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(norm = norm),
                                                   rowData = S4Vectors::DataFrame(rnames),
                                                   colData = S4Vectors::DataFrame(cnames))
  dt <- convert_se_assay_to_dt(se = se, assay_name = "norm", include_metadata = FALSE)
  merged <- base::merge(df, S4Vectors::DataFrame(dt[, c("rId", "cId", "values")]))
  
  expect_equal(merged$r, merged$rId)
  expect_equal(merged$c, merged$cId)
  
  dt <- convert_se_assay_to_dt(se = se, assay_name = "norm", include_metadata = TRUE, retain_nested_rownames = TRUE)
  
  # Properly handles nested rownames.
  expect_equal(rownames(dt), as.character(nested_rnames))
  expect_true("norm_rownames" %in% colnames(dt))
  
  merged <- base::merge(df, S4Vectors::DataFrame(dt[, c("rnames", "cnames", "values")]))
  expect_equal(merged$r, merged$rnames)
  expect_equal(merged$c, merged$cnames)
})


test_that("merge_metrics argument of assay_to_dt works as expected", {
  headers <- get_header("response_metrics")
  m <- 20
  n <- length(headers)
  metrics <- as.data.frame(matrix(runif(m * 2 * n), nrow = m * 2, ncol = n))
  colnames(metrics) <- headers
  metrics$dr_metric <- rep(c("RV", "GR"), m)
  rId <- as.factor(rep(rep(seq(4), each = 2), 5))
  cId <- as.factor(rep(rep(seq(5), each = 2), 4))
  mat <- BumpyMatrix::splitAsBumpyMatrix(metrics, row = rId, column = cId)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(Metrics = mat),
                                                   rowData = S4Vectors::DataFrame(rownames(mat)),
                                                   colData = S4Vectors::DataFrame(colnames(mat)))
  
  obs <- convert_se_assay_to_dt(se, "Metrics")

  expect_equal(nrow(obs), m * 2)
  expect_true(all(c(colnames(get_header("metrics_names"))) %in% colnames(obs)))
  
  # Insert random column. 
  metrics2 <- metrics
  extra_col <- "SERENA_WILLIAMS"
  extra_val <- rep_len(LETTERS, nrow(metrics2))
  metrics2[[extra_col]] <- extra_val
  mat2 <- BumpyMatrix::splitAsBumpyMatrix(metrics2, row = rId, column = cId)
  se2 <- SummarizedExperiment::SummarizedExperiment(assays = list(Metrics = mat2),
                                                    rowData = S4Vectors::DataFrame(rownames(mat2)),
                                                    colData = S4Vectors::DataFrame(colnames(mat2)))
  
  obs2 <- convert_se_assay_to_dt(se2, "Metrics")
  expect_true(extra_col %in% colnames(obs2))
  expect_equal(metrics2[[extra_col]], extra_val)
  expect_true(all(c(colnames(get_header("metrics_names"))) %in% colnames(obs2)))
})
