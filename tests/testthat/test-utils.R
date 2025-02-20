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
  expect_warning(out <- .clean_key_inputs(keys, cols))
  expect_equal(out, cols)
})


test_that("assert_equal_input_len works as expected", {
  ec50 <- 0.5
  x_0 <- 1
  x_inf <- 0.1
  h <- 2
  efficacy <- 0.6
  expect_equal(assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h), NULL)
  
  # Bad lengths.
  ec50 <- c(0.5, 0.5)
  expect_error(assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h))
  
  # Length 1 fit parameters.
  ec50 <- 0.5
  efficacy <- c(0.6, 0.7, 0.8)
  expect_equal(assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h), NULL)
  
  # Length 1 outlier.
  ec50 <- c(0.5, 0.6)
  x_0 <- c(1, 0.9)
  x_inf <- c(0.1, 0.15)
  h <- c(2, 2)
  efficacy <- 0.6
  expect_equal(assert_equal_input_len(outlier = efficacy, ec50, x_0, x_inf, h), NULL)
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
  expect_length(v1, 9)
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
  ligand_data <- get_synthetic_data("finalMAE_wLigand")
  metrics_data <- convert_se_assay_to_dt(ligand_data[[1]], "Metrics")
  data.table::setnames(metrics_data,
                       prettify_flat_metrics(names(metrics_data),
                                             human_readable = TRUE))
  avg_metrics_data <- average_biological_replicates_dt(dt = metrics_data,
                                                       var = "Ligand",
                                                       prettified = TRUE)
  expect_equal(dim(metrics_data), c(60, 29))
  expect_equal(dim(avg_metrics_data), c(40, 28))
  expect_true(!"Ligand" %in% names(avg_metrics_data))
  
  avg_metrics_data2 <- average_biological_replicates_dt(dt = metrics_data,
                                                        var = "Ligand",
                                                        prettified = TRUE,
                                                        add_sd = TRUE)
  
  expect_equal(dim(avg_metrics_data2), c(40, 44))
  expect_equal(sum(grepl("_sd", names(avg_metrics_data2))), 15)
  expect_true("count" %in% names(avg_metrics_data2))
  
  # protection against regression
  # fit_type correctly recognized in wide and long format
  sdata <- get_synthetic_data("finalMAE_small")
  smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
  tdata <- smetrics_data[1:8, ]
  tdata$Gnumber <- tdata$Gnumber[1]
  tdata$DrugName <- tdata$DrugName[1]
  tdata$source_id <- paste0("DS", rep(1:4, each = 2))
  tdata$fit_type <- letters[1:8]
  
  av1b <- average_biological_replicates_dt(tdata, var = "source_id")
  av1f <- flatten(
    av1b,
    groups = c("normalization_type", "fit_source"),
    wide_cols = get_header("response_metrics")
  )
  
  av2f <- flatten(
    tdata,
    groups = c("normalization_type", "fit_source"),
    wide_cols = get_header("response_metrics")
  )
  av2b <- average_biological_replicates_dt(av2f, var = "source_id")
  expect_true(all.equal(av1f, av2b))
  expect_true(nrow(av1f) == 1)
  av1i <- average_biological_replicates_dt(tdata, var = "source_id", fit_type_average_fields = "bad_value")
  expect_true(nrow(av1i) == 8)
  
})

test_that("has_single_codrug_data works as expected", {
  expect_false(has_single_codrug_data("un_col"))
  expect_true(has_single_codrug_data(get_prettified_identifiers(c(
    "concentration2", "drug_name2"
  ), simplify = FALSE)))
  expect_true(
    has_single_codrug_data(c("Concentration 2", "Drug Name 2", "anything")))
  expect_true(
    has_single_codrug_data(c("Concentration_2", "DrugName_2", "anything"), 
                           prettify_identifiers = FALSE))
  expect_true(
    has_single_codrug_data(c("Concentration 3", "Drug Name 3", "tissue"),
                           codrug_identifiers = c("concentration3", "drug_name3")))
  expect_true(
    has_single_codrug_data(c("Concentration_3", "DrugName_3", "tissue"),
                           prettify_identifiers = FALSE,
                           codrug_identifiers = c("concentration3", "drug_name3")))
  
  expect_error(
    has_single_codrug_data(list(drug = "test")),
    "Assertion on 'cols' failed: Must be of type 'character', not 'list'."
  )
  expect_error(
    has_single_codrug_data(c("Concentration 2", "Drug Name 2"), prettify_identifiers = "str"),
    "Assertion on 'prettify_identifiers' failed: Must be of type 'logical flag', not 'character'."
  )
  expect_error(
    has_single_codrug_data(c("drug", "conc"), codrug_identifiers = c(1, 2)),
    "Assertion on 'all(codrug_identifiers %in% names(get_env_identifiers(simplify = TRUE)))'",
    fixed = TRUE
  )
})

test_that("has_valid_codrug_data works as expected", {
  dt1 <-
    data.table::data.table(
      "Drug Name" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "Drug Name 2" = "untreated",
      "Concentration 2" = NA,
      "Drug Name 3" = "untreated",
      "Concentration 3" = NA
    )
  dt2 <-
    data.table::data.table(
      "Drug Name" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "Drug Name 2" = letters[4:6],
      "Concentration 2" = 4:6,
      "Drug Name 3" = letters[7:9],
      "Concentration 3" = 7:9
    )
  
  dt3 <-
    data.table::data.table(
      "DrugName" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "DrugName_2" = letters[4:6],
      "Concentration_2" = 4:6,
      "DrugName_3" = letters[7:9],
      "Concentration_3" = 7:9
    )
  
  expect_true(has_valid_codrug_data(dt2))
  expect_false(has_valid_codrug_data(dt2, prettify_identifiers = FALSE))
  expect_true(
    has_valid_codrug_data(
      dt2,
      codrug_name_identifier = "drug_name3",
      codrug_conc_identifier = "concentration3"
    )
  )
  
  expect_false(has_valid_codrug_data(dt3))
  expect_true(has_valid_codrug_data(dt3, prettify_identifiers = FALSE))
  expect_true(
    has_valid_codrug_data(
      dt3,
      prettify_identifiers = FALSE,
      codrug_name_identifier = "drug_name3",
      codrug_conc_identifier = "concentration3"
    )
  )
  
  expect_false(has_valid_codrug_data(dt2[, c("Drug Name", "Concentration")]))
  
  dt2[["Concentration 2"]] <- NA
  expect_false(has_valid_codrug_data(dt2))
  
  dt2[["Drug Name 3"]] <- "untreated"
  expect_false(
    has_valid_codrug_data(
      dt2,
      codrug_name_identifier = "drug_name3",
      codrug_conc_identifier = "concentration3"
    )
  )
  
  expect_error(
    has_valid_codrug_data(colnames(dt1)),
    "Assertion on 'data' failed: Must be a data.table, not character."
  )
  expect_error(
    has_valid_codrug_data(dt1, prettify_identifiers = "str"),
    "Assertion on 'prettify_identifiers' failed: Must be of type 'logical flag', not 'character'."
  )
  expect_error(
    has_valid_codrug_data(dt1, codrug_name_identifier = c("id1", "id2")),
    "Assertion on 'codrug_name_identifier' failed: Must have length 1."
  )
  expect_error(
    has_valid_codrug_data(dt1, codrug_conc_identifier = c("id1", "id2")),
    "Assertion on 'codrug_conc_identifier' failed: Must have length 1."
  )
  
})

test_that("remove_codrug_data works as expected", {
  dt1 <-
    data.table::data.table(
      "Drug Name" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "Drug Name 2" = "untreated",
      "Concentration 2" = NA,
      "Drug Name 3" = "untreated",
      "Concentration 3" = NA
    )
  
  sdt <- remove_codrug_data(dt1)
  exp_cols <- c("Drug Name", "Concentration", "Drug Name 3", "Concentration 3")
  expect_identical(colnames(sdt), exp_cols)
  
  sdt <- remove_codrug_data(dt1, codrug_identifiers = c("drug_name3", "concentration3"))
  exp_cols <- c("Drug Name", "Concentration", "Drug Name 2", "Concentration 2")
  expect_identical(colnames(sdt), exp_cols)
  
  dt2 <-
    data.table::data.table(
      "DrugName" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "DrugName_2" = "untreated",
      "Concentration_2" = NA,
      "DrugName_3" = "untreated",
      "Concentration_3" = NA
    )
  
  sdt <- remove_codrug_data(dt2, prettify_identifiers = FALSE)
  exp_cols <-  c("DrugName", "Concentration", "DrugName_3", "Concentration_3")
  expect_identical(colnames(sdt), exp_cols)
  
  expect_error(
    remove_codrug_data(colnames(dt1)),
    "Assertion on 'data' failed: Must be a data.table, not character."
  )
  expect_error(
    remove_codrug_data(dt1, prettify_identifiers = 1),
    "Assertion on 'prettify_identifiers' failed: Must be of type 'logical flag', not 'double'."
  )
  expect_error(
    remove_codrug_data(dt1, codrug_identifiers = "str"),
    "failed: Must be TRUE."
  )
})

test_that("is_combo_data works fine", {
  rdata <- data.table::data.table(Gnumber = seq_len(10),
                                  Concentration = runif(10), Ligand = c(rep(0.5, 5), rep(0, 5)))
  se <- SummarizedExperiment::SummarizedExperiment(rowData = rdata)
  expect_false(is_combo_data(se))
  
  nrows <- 10
  ncols <- 6
  mx <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  se <- SummarizedExperiment::SummarizedExperiment(
    rowData = rdata, 
    assays = list(excess = mx, scores = mx, isobolograms = mx))
  expect_true(is_combo_data(se))
  
  expect_error(is_combo_data(list()), "Must inherit from class 'SummarizedExperiment'")
})

test_that("get_additional_variables works as expected", {
  rdata1 <- data.table::data.table(Gnumber = seq_len(10),
                                   Concentration = runif(10), 
                                   Ligand = c(rep(0.5, 5), rep(0, 5)))
  rdata2 <- data.table::data.table(`Drug Name` = seq_len(10),
                                   Concentration = runif(10), 
                                   Ligand = c(rep(0.5, 10)),
                                   Replicate = seq_len(2))
  rdata3 <- data.table::data.table(Gnumber = seq_len(10),
                                   Concentration = runif(10))
  rdata4 <- data.table::data.table(Gnumber = seq_len(10),
                                   Concentration = runif(10), 
                                   Concentration_2 = runif(10), 
                                   Ligand = c(rep(0.5, 10)),
                                   Ligand = c(rep(0.1, 5), rep(0, 5)))
  
  add_var1 <- get_additional_variables(rdata1, prettified = TRUE)
  add_var2_nonunique <- get_additional_variables(rdata2, prettified = TRUE)
  add_var2_unique <- get_additional_variables(rdata2, unique = TRUE, prettified = TRUE)
  add_var3 <- get_additional_variables(rdata3, prettified = TRUE)
  add_var4_nonunique <- get_additional_variables(rdata4, prettified = TRUE)
  add_var4_unique <- get_additional_variables(rdata4, unique = TRUE, prettified = TRUE)
  
  expect_equal(add_var1, "Ligand")
  expect_equal(add_var2_nonunique, "Replicate")
  expect_equal(add_var2_unique, c("Ligand", "Replicate"))
  expect_equal(add_var3, NULL)
  expect_equal(add_var4_nonunique, "Concentration_2")
  expect_equal(add_var4_unique, c("Concentration_2", "Ligand"))
  
  expect_equal(get_additional_variables(FALSE), NULL)
  expect_equal(get_additional_variables(NA), NULL)
  expect_equal(get_additional_variables(c(1, 2, 3)), NULL)
  expect_equal(get_additional_variables(unlist(rdata1)), NULL)
  
  
  rdata5 <- data.table::data.table(Gnumber = seq_len(10),
                                   Concentration = runif(10), 
                                   Concentration_2 = runif(10), 
                                   `IC50 (GDS)` = runif(10))
  expect_equal(get_additional_variables(rdata5), NULL)
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

test_that("calc_sd works as expected", {
  expect_equal(calc_sd(c(1, 2, 3, 4, 5)), sd(c(1, 2, 3, 4, 5), na.rm = TRUE))
  expect_equal(calc_sd(c(10, 20, 30)), sd(c(10, 20, 30), na.rm = TRUE))
  expect_equal(calc_sd(c(1)), 0)
  expect_true(is.na(calc_sd("2")))
  expect_true(is.na(calc_sd(TRUE)))
  expect_true(is.na(calc_sd(numeric(0))))
  expect_equal(calc_sd(c(1, 2, NA, 4, 5)), sd(c(1, 2, NA, 4, 5), na.rm = TRUE))
  expect_true(is.na(calc_sd(c(NA, NA, NA))))
})

test_that("remove_drug_batch", {
  # no suffix - nothing changes
  expect_equal(remove_drug_batch("G00060245"), "G00060245")
  # expected suffix - remove
  expect_equal(remove_drug_batch("G00060245.1-8"), "G00060245")
  expect_equal(remove_drug_batch("G00060245.18"), "G00060245")
  expect_equal(remove_drug_batch("G02948263.1-1.DMA"), "G02948263")
  # (single codrug) - remove
  expect_equal(remove_drug_batch("G03252046.1-2;G00376771"), "G03252046")
  # (two codrugs) - remove
  expect_equal(
    remove_drug_batch("G03256376.1-2;G00376771.1-19;G02557755"), "G03256376")
  
  # (Gnumber followed by the ",") -remove
  expect_equal(remove_drug_batch("G00018838, Cisplatin"), "G00018838")
  
  # suffix added by set_unique_drug_names_dt function (prevent duplication) - nothing changes
  expect_equal(remove_drug_batch("G00060245_(G00060245.1-8)"),
               "G00060245_(G00060245.1-8)")
  
  # test non-default values of other parameters
  expect_equal(remove_drug_batch("DRUG_01.123", drug_p = "DRUG_[0-9]+"),
               "DRUG_01")
  expect_equal(remove_drug_batch("G00001234:22-1", sep_p = ":"), "G00001234")
  expect_equal(remove_drug_batch("G00001234.28", batch_p = "[0-9]+"),
               "G00001234")
  
  expect_error(remove_drug_batch(list(drug = "G00000001")), "Assertion on 'drug_vec' failed")
  expect_error(remove_drug_batch("G00000001", drug_p = list(1)),
               "Assertion on 'drug_p' failed")
  expect_error(remove_drug_batch("G00000001", sep_p = list(1)),
               "Assertion on 'sep_p' failed")
  expect_error(remove_drug_batch("G00000001", batch_p = list(1)),
               "Assertion on 'batch_p' failed")
  
})

test_that("cap_assay_infinities", {
  # single-agent data - data expected tests
  sdata <- get_synthetic_data("finalMAE_small")
  smetrics_data <- convert_se_assay_to_dt(sdata[[get_supported_experiments("sa")]], 
                                          "Metrics")
  ## add some Infs/-Infs
  smetrics_data$xc50[1:3] <- -Inf
  smetrics_data$xc50[100:103] <- Inf
  saveraged_data <- convert_se_assay_to_dt(sdata[[get_supported_experiments("sa")]], 
                                           "Averaged")
  smetrics_data2 <- cap_assay_infinities(saveraged_data, 
                                         smetrics_data, 
                                         experiment_name = get_supported_experiments("sa")) # default
  smetrics_data3 <- cap_assay_infinities(saveraged_data,
                                         smetrics_data,
                                         experiment_name = get_supported_experiments("sa"),
                                         capping_fold = 1)
  
  ## data with inf/-inf values
  inf_idx <- which(is.infinite(smetrics_data$xc50))
  expect_true(NROW(inf_idx) > 0)
  # no Inf/-Inf after running the function
  inf_idx2 <- which(is.infinite(smetrics_data2$xc50))
  expect_true(NROW(inf_idx2) == 0)
  ##  Inf values
  expect_identical(max(smetrics_data2[inf_idx, "xc50"]), 5 * max(smetrics_data3[inf_idx, "xc50"]))
  expect_identical(min(smetrics_data2[inf_idx, "xc50"]), min(smetrics_data3[inf_idx, "xc50"]) / 5)
  
  ## data without infinities
  smetrics_data4 <- cap_assay_infinities(saveraged_data, 
                                         smetrics_data2, 
                                         experiment_name = get_supported_experiments("sa"))
  expect_identical(smetrics_data2, smetrics_data4)
  
  ## non-default column to be changed
  smetrics_data5 <- smetrics_data
  smetrics_data5$custom_col <- smetrics_data5$xc50
  smetrics_data6 <- cap_assay_infinities(saveraged_data,
                                         smetrics_data5,
                                         experiment_name = "single-agent",
                                         col = "custom_col")
  expect_identical(smetrics_data2$xc50, smetrics_data6$custom_col)
  expect_true(any(smetrics_data6$xc50 != smetrics_data2$xc50))
  
  # combination data - data expected tests
  cdata <- get_synthetic_data("finalMAE_combo_matrix_small")
  scaveraged_data <- convert_se_assay_to_dt(cdata[[get_supported_experiments("combo")]], 
                                            "Averaged")
  scmetrics_data <- convert_se_assay_to_dt(cdata[[get_supported_experiments("combo")]], 
                                           "Metrics")
  scmetrics_data2 <- cap_assay_infinities(scaveraged_data, 
                                          scmetrics_data, 
                                          experiment_name = get_supported_experiments("combo"))
  scmetrics_data3 <- cap_assay_infinities(scaveraged_data,
                                          scmetrics_data,
                                          experiment_name = get_supported_experiments("combo"),
                                          capping_fold = 1)
  
  ## data with inf/-inf values
  infc_idx <- which(is.infinite(scmetrics_data$xc50))
  expect_true(NROW(infc_idx) > 0)
  ## no Inf/-Inf after running the function
  infc_idx2 <- which(is.infinite(scmetrics_data2$xc50))
  expect_true(NROW(infc_idx2) == 0)
  ##  Inf values
  expect_identical(max(scmetrics_data2[infc_idx, "xc50"]), 5 * max(scmetrics_data3[infc_idx, "xc50"]))
  expect_identical(min(scmetrics_data2[infc_idx, "xc50"]), min(scmetrics_data3[infc_idx, "xc50"]) / 5)
  
  ## data without infinities
  scmetrics_data4 <- cap_assay_infinities(scaveraged_data, 
                                          scmetrics_data2, 
                                          experiment_name = get_supported_experiments("combo"))
  expect_identical(scmetrics_data2, scmetrics_data4)
  
  # test non-default values of other parameters
  expect_error(cap_assay_infinities(list(a = 2)), "Must be a data.table")
  expect_error(cap_assay_infinities(saveraged_data, list(a = 2)),
               "Must be a data.table")
  expect_error(cap_assay_infinities(saveraged_data, 
                                    smetrics_data, 
                                    experiment_name = "test"),
               "Must be element of set ")
  expect_error(cap_assay_infinities(saveraged_data,
                                    smetrics_data,
                                    experiment_name = get_supported_experiments("sa"),
                                    col = 2),
               "Must be of type 'string'")
  expect_error(cap_assay_infinities(saveraged_data,
                                    smetrics_data,
                                    experiment_name = get_supported_experiments("sa"),
                                    col = "no_col"),
               "Must be element of set")
  expect_error(cap_assay_infinities(saveraged_data,
                                    smetrics_data,
                                    experiment_name = get_supported_experiments("sa"),
                                    capping_fold = "x"),
               "Must be of type 'number'")
})
