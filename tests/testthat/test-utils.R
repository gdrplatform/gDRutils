### input MAEs ###
empty_se <-
  SummarizedExperiment::SummarizedExperiment(assays = list(data.frame()))
empty_mae <-
  MultiAssayExperiment::MultiAssayExperiment(experiments = MultiAssayExperiment::ExperimentList(test = empty_se))
maeReal <-
  get_synthetic_data("finalMAE_combo_2dose_nonoise2")

partially_empty_mae <-
  MultiAssayExperiment::MultiAssayExperiment(experiments = (MultiAssayExperiment::ExperimentList(
    experiments = c(
      MultiAssayExperiment::experiments(empty_mae),
      MultiAssayExperiment::experiments(maeReal)
    )
  )))

test_that(".clean_key_inputs works as expected", {
  keys <- LETTERS[1:5]
  cols <- LETTERS[1:3]
  expect_warning(out <- gDRutils:::.clean_key_inputs(keys, cols))
  expect_equal(out, cols)
})


test_that("assert_equal_input_len works as expected", {
  ec50 <- 0.5
  x_0 <- 1
  x_inf <- 0.1
  h <- 2
  efficacy <- 0.6
  expect_equal(gDRutils:::assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h), NULL)

  # Bad lengths.
  ec50 <- c(0.5, 0.5)
  expect_error(gDRutils:::assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h))

  # Length 1 fit parameters.
  ec50 <- 0.5
  efficacy <- c(0.6, 0.7, 0.8)
  expect_equal(gDRutils:::assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h), NULL)

  # Length 1 outlier.
  ec50 <- c(0.5, 0.6)
  x_0 <- c(1, 0.9)
  x_inf <- c(0.1, 0.15)
  h <- c(2, 2)
  efficacy <- 0.6
  expect_equal(gDRutils:::assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h), NULL)
})

test_that("assert_choices",  {
  ### expected values
  expect_null(assert_choices(letters[1], letters))
  expect_null(assert_choices(letters[1:2], letters))
  expect_null(assert_choices(1:5, 1:10))
  expect_null(assert_choices(1, 1:10))

  ### errors
  err_msg <- sprintf("Assertion on '%s' failed.", letters[8])
  expect_error(assert_choices(letters[5:8], letters[1:7]), err_msg)
  err_msg <-
    sprintf("Assertion on '%s, %s' failed.", letters[7], letters[8])
  expect_error(assert_choices(letters[5:8], letters[1:6]), err_msg)
})

test_that("MAEpply works as expected", {
  list1 <- MAEpply(maeReal, SummarizedExperiment::assayNames)
  expect_length(list1, 2)
  expect_true(inherits(list1, "list"))
  v1 <- unique(MAEpply(maeReal, SummarizedExperiment::assayNames, unify = TRUE))
  expect_length(v1, 14)
  expect_true(inherits(v1, "character"))

  v2 <- unique(MAEpply(maeReal, SummarizedExperiment::rowData, unify = TRUE))
  expect_equal(dim(v2), c(7, 7))
  checkmate::expect_data_table(v2)
})

test_that("is_mae_empty works as expected", {
  expect_false(is_mae_empty(maeReal))
  expect_true(is_mae_empty(empty_mae))
})

test_that("is_any_exp_empty works as expected", {
  expect_false(is_any_exp_empty(maeReal))
  expect_true(is_any_exp_empty(empty_mae))
  expect_true(is_any_exp_empty(partially_empty_mae))
})

test_that("is_exp_empty works as expected", {
  expect_false(is_exp_empty(maeReal[[1]]))
  expect_true(is_exp_empty(empty_mae[[1]]))
})

test_that("get_non_empty_assays works as expected", {
  expect_identical(get_non_empty_assays(maeReal), get_non_empty_assays(partially_empty_mae))
  expect_identical(get_non_empty_assays(empty_mae), character(0))
})

test_that("mrowData works as expected", {
  mr <- mrowData(maeReal)
  expect_equal(dim(mr), c(7, 7))
  checkmate::expect_data_table(mr)

  mr <- mrowData(empty_mae)
  expect_identical(mr, data.table::data.table())
})

test_that("mcolData works as expected", {
  mc <- mcolData(maeReal)
  expect_equal(dim(mc), c(6, 4))
  checkmate::expect_data_table(mc)

  mc <- mcolData(empty_mae)
  expect_identical(mc, data.table::data.table())
})

test_that("apply_bumpy_function works as expected", {
  n <- 1000
  df <- data.table::data.table(row = sample(LETTERS, n, replace = TRUE),
    column = sample(LETTERS, n, replace = TRUE),
    a = runif(n),
    b = runif(n))
  bumpy <- BumpyMatrix::splitAsBumpyMatrix(df[, c("a", "b")], row = df$row, column = df$column)
  se <- SummarizedExperiment::SummarizedExperiment(assays =
    list(bumpy = bumpy))

  # Assertions.
  expect_error(apply_bumpy_function(se, req_assay_name = "nonexistent", out_assay_name = "misc"),
    regex = "'nonexistent' is not on of the available assays: 'bumpy'")

  # Output is bumpy matrix.
  FUN <- function(x) {
    data.table::data.table(y = x$a + x$b, z = x$a - x$b)
  }
  bumpy_out <- apply_bumpy_function(se, FUN = FUN, req_assay_name = "bumpy", out_assay_name = "bumpy_mtx")
  expect_true(is(bumpy_out, "SummarizedExperiment"))
  expect_true("bumpy_mtx" %in% SummarizedExperiment::assayNames(bumpy_out))

  bumpy_in_df <- BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(se, "bumpy"))
  bumpy_in_df$y <- bumpy_in_df$a + bumpy_in_df$b
  bumpy_in_df$z <- bumpy_in_df$a - bumpy_in_df$b

  bumpy_out_df <- BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(bumpy_out, "bumpy_mtx"))
  keep_cols <- c("row", "column", "y", "z")
  expect_equal(sort(bumpy_in_df[, keep_cols]), sort(bumpy_out_df[, keep_cols]))
})

test_that("loop works as expected", {
  n <- 10
  listRunif <- lapply(seq_len(n), runif)
  sumOfList <- loop(listRunif, sum)
  expect_true(is(sumOfList, "list"))
  expect_length(unlist(sumOfList), n)
})

test_that("get_synthetic_data works as expected", {
  expect_is(get_synthetic_data("finalMAE_small.qs"),
            "MultiAssayExperiment")
  expect_is(get_synthetic_data("finalMAE_small"),
            "MultiAssayExperiment")
  expect_is(get_synthetic_data("finalMAE_small.RDS"),
            "MultiAssayExperiment")
  expect_error(get_synthetic_data("finalMAE_small.wrong"), "Failed to open")
})


test_that("geometric_mean works as expected", {
  expect_error(
    geometric_mean(x = "NULL"),
    "Assertion on 'x' failed: Must be of type 'numeric', not 'character'."
  )
  expect_error(
    geometric_mean(x = 1, fixed = "NULL"),
    "Assertion on 'fixed' failed: Must be of type 'logical flag', not 'character'."
  )
  expect_error(
    geometric_mean(x = 1, maxlog10Concentration = "NULL"),
    "Assertion on 'maxlog10Concentration' failed: Must be of type 'numeric', not 'character'."
  )
  
  expect_equal(geometric_mean(c(2, 8)), 4)
  expect_equal(geometric_mean(c(0.02, 8)), 0.4)
  
  expect_equal(round(geometric_mean(c(0.000000000002, 8)), digits = 5), 0.00894)
  expect_equal(round(geometric_mean(c(0.000000000002, 8), fixed = TRUE), digits = 5), 0.00894)
  expect_equal(geometric_mean(c(0.000000000002, 8), fixed = FALSE), 0.000004)
  
  expect_equal(geometric_mean(c(2, 800)), 10)
  expect_equal(geometric_mean(c(2, 800), fixed = TRUE), 10)
  expect_equal(geometric_mean(c(2, 800), fixed = FALSE), 40)
  
  expect_equal(round(
    geometric_mean(c(2, 8), fixed = TRUE, maxlog10Concentration = 1),
    digits = 5
  ), 4)
  expect_equal(round(
    geometric_mean(c(2, 8), fixed = TRUE, maxlog10Concentration = 0.1),
    digits = 5
  ), 3.54813)
  
})


test_that("average_biological_replicates_dt works as expected", {
  ligand_data <- gDRutils::get_synthetic_data("finalMAE_wLigand")
  metrics_data <- convert_se_assay_to_dt(ligand_data[[1]], "Metrics")
  data.table::setnames(metrics_data,
                       prettify_flat_metrics(names(metrics_data),
                                             human_readable = TRUE))
  avg_metrics_data <- average_biological_replicates_dt(dt = metrics_data,
                                                       var = "Ligand")
  expect_equal(dim(metrics_data), c(60, 27))
  expect_equal(dim(avg_metrics_data), c(40, 26))
  expect_true(!"Ligand" %in% names(avg_metrics_data))
})

test_that("get_duplicated_rows works as expected", {
  
  # single column
  expect_equal(
    get_duplicated_rows(DF1co, col_names = "DrugName"),
    c(1, 2)
  )
  # single column with only duplicates
  expect_equal(
    get_duplicated_rows(DF1co, col_names = "DrugName_2"),
    c(1, 2, 3)
  )
  # single column without duplicates
  expect_equal(
    get_duplicated_rows(DF1co, col_names = c("Gnumber")),
    integer()
  )
  # multiple columns
  expect_equal(
    get_duplicated_rows(DF1co, col_names = c("DrugName_2", "DrugName")),
    c(1, 2)
  )
  
  expect_error(get_duplicated_rows(DF1co, c("DrugName", "Fake Column")),
               "Assertion on 'all(col_names %in% colnames(x))' failed: Must be TRUE.", fixed = TRUE)
})