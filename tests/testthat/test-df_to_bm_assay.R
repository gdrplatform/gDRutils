test_that("df_to_bm_assay works with default parameters", {
  dt <- data.table::data.table(Gnumber = "G1", clid = "CL1", x = 1)
  result <- df_to_bm_assay(dt)
  expect_s4_class(result, "BumpyMatrix")
})

test_that("df_to_bm_assay uses precomputed_metadata when provided", {
  dt <- data.table::data.table(Gnumber = "G1", clid = "CL1", x = 1)
  metadata <- split_SE_components(dt)
  result_default <- df_to_bm_assay(dt)
  result_precomputed <- df_to_bm_assay(dt, precomputed_metadata = metadata)
  expect_identical(dim(result_default), dim(result_precomputed))
  expect_identical(rownames(result_default), rownames(result_precomputed))
  expect_identical(colnames(result_default), colnames(result_precomputed))
})

test_that("df_to_bm_assay validates precomputed_metadata type", {
  dt <- data.table::data.table(Gnumber = "G1", clid = "CL1", x = 1)
  expect_error(
    df_to_bm_assay(dt, precomputed_metadata = "not_a_list"),
    "list"
  )
})

test_that("df_to_bm_assay falls back to recomputation for incomplete precomputed_metadata", {
  dt <- data.table::data.table(Gnumber = "G1", clid = "CL1", x = 1)
  incomplete <- list(condition_md = data.table::data.table(clid = "CL1"))
  result_fallback <- df_to_bm_assay(dt, precomputed_metadata = incomplete)
  result_default <- df_to_bm_assay(dt)
  expect_identical(dim(result_fallback), dim(result_default))
})

test_that("convert_dt_to_se builds an SE over all row and column ids", {
  dt <- data.table::CJ(compound = c("cA", "cB"), target = c("t1", "t2"),
                       Concentration = c(0.1, 1))
  dt[, Activity := seq_len(.N)]
  dt <- dt[!(compound == "cB" & target == "t2")]

  se <- convert_dt_to_se(dt, "compound", "target", "Biochemical")

  expect_identical(dimnames(se), list(c("cA", "cB"), c("t1", "t2")))
  expect_identical(SummarizedExperiment::assayNames(se), "Biochemical")
  expect_identical(SummarizedExperiment::rowData(se)$compound, c("cA", "cB"))
  expect_identical(SummarizedExperiment::colData(se)$target, c("t1", "t2"))
  mx <- SummarizedExperiment::assay(se, "Biochemical")
  expect_identical(colnames(mx["cA", "t1"][[1]]), c("Concentration", "Activity"))
  expect_equal(mx["cA", "t2"][[1]]$Activity, dt[compound == "cA" & target == "t2", Activity])
  expect_equal(NROW(mx["cB", "t2"][[1]]), 0L)

  # several id columns are pasted into one name and kept apart in rowData
  dt[, batch := "b1"]
  se2 <- convert_dt_to_se(dt, c("compound", "batch"), "target", "Biochemical")
  expect_identical(rownames(se2), c("cA_b1", "cB_b1"))
  expect_identical(SummarizedExperiment::rowData(se2)$batch, c("b1", "b1"))
})

test_that("convert_dt_to_se rejects ambiguous or empty input", {
  dt <- data.table::data.table(a = c("x_y", "x"), b = c("z", "y_z"), target = "t1", v = 1:2)
  expect_error(convert_dt_to_se(dt, c("a", "b"), "target", "B"), "paste to the same name: x_y_z")
  expect_error(convert_dt_to_se(dt, "a", "a", "B"))
  expect_error(convert_dt_to_se(dt[, .(a, target)], "a", "target", "B"), "no columns besides")
})

test_that("convert_dt_to_se output goes through apply_fit", {
  dt <- data.table::CJ(compound = c("cA", "cB"), target = "t1",
                       Concentration = c(0.1, 1), assay_type = "IC50")
  dt[, Activity := c(90, 40, 80, 20)]
  se <- convert_dt_to_se(dt, "compound", "target", "Biochemical")
  register_fit_profile("test_convert_dt", slicing_cols = "assay_type", input_assay = "Biochemical")
  on.exit(rm("test_convert_dt", envir = .fit_profile_env), add = TRUE)

  out <- apply_fit(se, function(dt) list(mean_activity = mean(dt$Activity)),
                   "test_convert_dt", output_assay = "Fits", fit_source = "test")
  res <- convert_se_assay_to_dt(out, "Fits")
  expect_equal(res[order(compound), mean_activity], c(65, 50))
})
