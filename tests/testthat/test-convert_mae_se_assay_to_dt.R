test_that("convert_se_assay_to_dt works as expected", {
  m <- 20
  n <- 10
  rnames <- LETTERS[1:m]
  cnames <- letters[1:n]
  
  # Normal matrix.
  ref_gr_value <- matrix(runif(m * n), nrow = m, ncol = n, dimnames = list(rnames, cnames))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(RefGRvalue = ref_gr_value),
                                                   rowData = S4Vectors::DataFrame(rnames),
                                                   colData = S4Vectors::DataFrame(cnames))
  
  dt <- convert_se_assay_to_dt(se = se, assay_name = "RefGRvalue", include_metadata = FALSE)
  expect_equal(dt$RefGRvalue, as.vector(ref_gr_value))
  expect_equal(dim(dt), c(m * n, 3))
  
  dt <- convert_se_assay_to_dt(se = se, assay_name = "RefGRvalue", include_metadata = TRUE)
  
  expect_equal(data.table::setorder(dt, cId)[["RefGRvalue"]], as.vector(ref_gr_value))
  expect_equal(dim(dt), c(m * n, 5))
  expect_equal(dt$rnames, as.character(dt$rId))
  expect_equal(dt$cnames, as.character(dt$cId))
  
  # BumpyDataFrameMatrix.
  df <- S4Vectors::DataFrame(r = rep(rnames, n), c = rep(cnames, m), values = runif(m * n))
  nested_rnames <- paste0("A", seq(nrow(df)))
  rownames(df) <- nested_rnames
  
  norm <- BumpyMatrix::splitAsBumpyMatrix(df, row = df$r, column = df$c)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(norm = norm),
                                                   rowData = S4Vectors::DataFrame(rnames),
                                                   colData = S4Vectors::DataFrame(cnames))
  dt <- convert_se_assay_to_dt(se = se, assay_name = "norm", include_metadata = FALSE)
  merged <- base::merge(df, S4Vectors::DataFrame(dt[, c("rId", "cId", "values")]))
  
  expect_equal(merged$r, merged$rId)
  expect_equal(merged$c, merged$cId)
  
  dt <- convert_se_assay_to_dt(se = se, assay_name = "norm", 
                               include_metadata = TRUE, retain_nested_rownames = TRUE)
  
  # Properly handles nested rownames.
  expect_equal(sort(dt$norm_rownames), sort(as.character(nested_rnames)))
  expect_true("norm_rownames" %in% colnames(dt))
  
  merged <- base::merge(df, S4Vectors::DataFrame(dt[, c("rnames", "cnames", "values")]))
  expect_equal(merged$r, merged$rnames)
  expect_equal(merged$c, merged$cnames)
})


test_that("merge_metrics argument of assay_to_dt works as expected", {
  headers <- get_header("response_metrics")
  m <- 20
  n <- length(headers)
  metrics <- data.table::as.data.table(matrix(runif(m * 2 * n), nrow = m * 2, ncol = n))
  data.table::setnames(metrics, headers)
  metrics[, dr_metric := rep(c("RV", "GR"), m)]
  rId <- as.factor(rep(rep(seq(4), each = 2), 5))
  cId <- as.factor(rep(rep(seq(5), each = 2), 4))
  mat <- BumpyMatrix::splitAsBumpyMatrix(metrics, row = rId, column = cId)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(Metrics = mat),
                                                   rowData = S4Vectors::DataFrame(rownames(mat)),
                                                   colData = S4Vectors::DataFrame(colnames(mat)))
  
  obs <- convert_se_assay_to_dt(se, "Metrics")
  
  expect_equal(nrow(obs), m * 2)
  expect_true(all(colnames(get_header("metrics_names")) %in% colnames(obs)))
  
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
  expect_true(all(colnames(get_header("metrics_names")) %in% colnames(obs2)))
})


test_that("convert_mae_assay_to_dt works as expected", {
  m <- 20
  n <- 10
  rnames <- LETTERS[1:m]
  cnames <- letters[1:n]
  
  # Normal matrix.
  ref_gr_value <-  matrix(runif(m * n), nrow = m, ncol = n, dimnames = list(rnames, cnames))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(RefGRvalue = ref_gr_value),
                                                   rowData = S4Vectors::DataFrame(rnames),
                                                   colData = S4Vectors::DataFrame(cnames))
  
  mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list("single-agent" = se))
  
  dt <- convert_mae_assay_to_dt(mae = mae, assay_name = "RefGRvalue", include_metadata = FALSE)
  checkmate::expect_data_table(dt)
  expect_equal(dt$RefGRvalue, as.vector(ref_gr_value))
  expect_equal(dim(dt), c(m * n, 3))
  
  dt <- suppressWarnings(convert_mae_assay_to_dt(mae = mae, assay_name = "RefGRvalue", 
                                                 include_metadata = FALSE, wide_structure = TRUE))
  checkmate::expect_data_table(dt)
  expect_equal(dt$RefGRvalue, as.vector(ref_gr_value))
  expect_equal(dim(dt), c(m * n, 3))
  expect_warning(convert_mae_assay_to_dt(mae = mae, assay_name = "RefGRvalue",
                                         include_metadata = FALSE, wide_structure = TRUE),
                 "assay is not class `BumpyMatrix`, wide_structure=TRUE ignored")
  
  dt <- convert_se_assay_to_dt(se = se, assay_name = "RefGRvalue", include_metadata = TRUE)
  expect_equal(data.table::setorder(dt, cId)[["RefGRvalue"]], as.vector(ref_gr_value))
  expect_equal(dim(dt), c(m * n, 5))
  expect_equal(dt$rnames, as.character(dt$rId))
  expect_equal(dt$cnames, as.character(dt$cId))
  
  se1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(RefGRvalue = ref_gr_value[1:10, ]),
    rowData = S4Vectors::DataFrame(rnames)[1:10, , drop = FALSE],
    colData = S4Vectors::DataFrame(cnames))
  
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(RefGRvalue = ref_gr_value[11:20, ]),
    rowData = S4Vectors::DataFrame(rnames)[11:20, , drop = FALSE],
    colData = S4Vectors::DataFrame(cnames))
  
  maeTwoExperiments <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list("single-agent" = se1,
                       "combination" = se2))
  
  dt1 <- convert_mae_assay_to_dt(mae = maeTwoExperiments, experiment_name = "single-agent",
                                 assay_name = "RefGRvalue", include_metadata = FALSE)
  checkmate::expect_data_table(dt1)
  expect_equal(dt1$RefGRvalue, as.vector(ref_gr_value[1:10, , drop = FALSE]))
  expect_equal(dim(dt1), c(m / 2 * n, 3))
  
  dt1 <- convert_mae_assay_to_dt(mae = maeTwoExperiments, experiment_name = "single-agent",
                                 assay_name = "RefGRvalue", include_metadata = TRUE)
  checkmate::expect_data_table(dt1)
  expect_equal(sort(data.table::setorder(dt1)[["RefGRvalue"]]), 
               sort(as.vector(ref_gr_value[1:10, ])))
  expect_equal(dim(dt1), c(m / 2 * n, 5))
  expect_equal(dt1$rnames, as.character(dt1$rId))
  expect_equal(dt1$cnames, as.character(dt1$cId))
  
  
  dt2 <- convert_mae_assay_to_dt(mae = maeTwoExperiments, experiment_name = "combination",
                                 assay_name = "RefGRvalue", include_metadata = FALSE)
  checkmate::expect_data_table(dt2)
  expect_equal(dt2$RefGRvalue, as.vector(ref_gr_value[11:20, , drop = FALSE]))
  expect_equal(dim(dt2), c(m / 2 * n, 3))
  
  dt2 <- convert_mae_assay_to_dt(mae = maeTwoExperiments, experiment_name = "combination",
                                 assay_name = "RefGRvalue", include_metadata = TRUE)
  checkmate::expect_data_table(dt2)
  expect_equal(sort(data.table::setorder(dt2)[["RefGRvalue"]]), 
               sort(as.vector(ref_gr_value[11:20, ])))
  expect_equal(dim(dt2), c(m / 2 * n, 5))
  expect_equal(dt2$rnames, as.character(dt2$rId))
  expect_equal(dt2$cnames, as.character(dt2$cId))
  
  
  dt3 <- convert_mae_assay_to_dt(mae = maeTwoExperiments,
                                 assay_name = "RefGRvalue", include_metadata = TRUE)
  checkmate::expect_data_table(dt3)
  expect_equal(data.table::setorder(dt3, cId)[["RefGRvalue"]], as.vector(ref_gr_value))
  expect_equal(dim(dt), c(m * n, 5))
  expect_equal(dt$rnames, as.character(dt$rId))
  expect_equal(dt$cnames, as.character(dt$cId))
  
  expect_warning(convert_mae_assay_to_dt(mae = maeTwoExperiments,
                                         assay_name = "Nonexistent"),
                 "assay 'Nonexistent' was not found in any of the following experiments")
  
  expect_warning(convert_mae_assay_to_dt(mae = maeTwoExperiments, assay_name = "RefGRvalue", 
                                         include_metadata = TRUE, wide_structure = TRUE),
                 "assay is not class `BumpyMatrix`, wide_structure=TRUE ignored")
  
  
  real_mae <- get_synthetic_data("finalMAE_small")
  
  dt_l <- convert_mae_assay_to_dt(mae = real_mae, assay_name = "Averaged", wide_structure = FALSE)
  checkmate::expect_data_table(dt_l)
  col_l <- c("normalization_type", "x", "x_std")
  expect_true(all(col_l %in% colnames(dt_l)))
  
  dt_w <- convert_mae_assay_to_dt(mae = real_mae, assay_name = "Averaged", wide_structure = TRUE)
  checkmate::expect_data_table(dt_w)
  col_w <- c("RelativeViability", "GRvalue", "std_RelativeViability", "std_GRvalue")
  expect_true(all(col_w %in% colnames(dt_w)))
  
  expect_equal(setdiff(colnames(dt_w), col_w), setdiff(colnames(dt_l), col_l))
  expect_equal(max(dt_l[, .N, by = "normalization_type"]$N), NROW(dt_w))
})


test_that("convert_se_assay_to_custom_dt works fine", {
  json_path <- system.file(package = "gDRutils", "test_settings_2.json")
  s <- get_settings_from_json(json_path = json_path)
  
  se <- get_synthetic_data("finalMAE_small")[[1]]
  dt1 <- convert_se_assay_to_custom_dt(se, assay_name = "Metrics")
  checkmate::expect_data_table(dt1, min.rows = 2, min.cols = 2)
  expect_true(all(s$METRIC_WISH_LIST %in% names(dt1)))
  dt2 <-
    convert_se_assay_to_custom_dt(se, assay_name = "Metrics", output_table = "Metrics_raw")
  checkmate::expect_data_table(dt2, min.rows = 2, min.cols = 2)
  expect_true(all(s$METRIC_WISH_LIST %in% names(dt2)))
  dt3 <-
    convert_se_assay_to_custom_dt(se, assay_name = "Metrics", output_table = "Metrics_initial")
  checkmate::expect_data_table(dt3, min.rows = 2, min.cols = 2)
  expect_false(identical(dt2, dt3))
  checkmate::expect_data_table(dt2, min.rows = 2, min.cols = 2)
  expect_true(all(c("x_mean", "x_AOC", "x_AOC_range", "xc50", "x_max", "ec50", 
                    "x_inf", "x_0", "h", "r2", "x_sd_avg", "fit_type") %in% names(dt3)))
  dt4 <- convert_se_assay_to_custom_dt(se, assay_name = "Averaged")
  checkmate::expect_data_table(dt2, min.rows = 2, min.cols = 2)
  expect_true(
    all(c("GR value", "Relative Viability", "Std GR value", "Std Relative Viability") %in% names(dt4)))
  dt5 <- convert_se_assay_to_custom_dt(se, assay_name = "Averaged", output_table = "Metrics")
  expect_true(identical(dt4, dt5))
  
  se2 <- get_synthetic_data("finalMAE_combo_matrix")[[1]]
  dt6 <- convert_se_assay_to_custom_dt(se2, assay_name = "Metrics")
  checkmate::expect_data_table(dt6, min.rows = 2, min.cols = 2)
  expect_true(all(s$METRIC_WISH_LIST %in% names(dt6)))
  dt7 <-
    convert_se_assay_to_custom_dt(se2, assay_name = get_combo_assay_names()[1])
  checkmate::expect_data_table(dt7, min.rows = 2, min.cols = 2)
  expect_true(all(names(get_combo_excess_field_names()) %in% names(dt7)))
  
  expect_error(convert_se_assay_to_custom_dt(as.list(se), assay_name = "Metrics"))
  expect_error(convert_se_assay_to_custom_dt(as.list(se), output_table = "Averaged"))
  expect_error(
    convert_se_assay_to_custom_dt(as.list(se), assay_name = "Averaged", output_table = "Metrics_raw"))
  expect_error(convert_se_assay_to_custom_dt(se, "xxx"))
  expect_error(convert_se_assay_to_custom_dt(se, "Metrics", "xxx"))
})



test_that("capVals works as expected", {
  dt1 <- data.table::data.table(
    `E Max` = c(-0.1, 0, 0.5, 1.2),
    `GR Max` = c(-1.1, -1, 0.5, 1.2),
    `RV AOC within set range` = c(-0.2, -0.1, 0, 3),
    `GR AOC within set range` = c(-0.2, -0.1, 0, 3),
    `GR50` = c(0, 1e-7, 10, 34),
    `IC50` = c(0, 1e-7, 10, 34),
    `EC50` = c(0, 1e-7, 10, 34),
    check.names = FALSE
  )
  dt2 <- data.table::data.table(
    `E Max` = c(0, 0, 0.5, 1.1),
    `GR Max` = c(-1, -1, 0.5, 1.1),
    `RV AOC within set range` = c(-0.1, -0.1, 0, 3),
    `GR AOC within set range` = c(-0.1, -0.1, 0, 3),
    `GR50` = c(1e-4, 1e-4, 10, 30),
    `IC50` = c(1e-4, 1e-4, 10, 30),
    `EC50` = c(NA, 1e-4, 10, 30),
    check.names = FALSE
  )
  dt3 <- data.table::data.table(
    A = LETTERS[1:10],
    B = letters[1:10],
    C = 1:10
  )
  dt1c <- capVals(dt1)
  dt2c <- capVals(dt2)
  dt2c_2 <- capVals(dt2[, 1:4])
  dt3c <- capVals(dt3)
  # remove index attribute, created inside capVals (for comparison purposes)
  attr(dt1c, "index") <- NULL
  attr(dt2c, "index") <- NULL
  attr(dt2c_2, "index") <- NULL
  attr(dt3c, "index") <- NULL
  
  expect_false(identical(dt1c, dt1))
  expect_identical(dt2c, dt2)
  expect_identical(dt2c_2, dt2[, 1:4])
  expect_identical(dt3c, dt3)
  
  # values are capped correctly
  expect_equal(dt1c, dt2)
  
  expect_error(capVals(as.list(dt1)), "Must be a data.table")
})

