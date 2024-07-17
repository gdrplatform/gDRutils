library(testthat)
context("combo-related functions")

test_that("convert_combo_data_to_dt", {
  test_mae <- get_synthetic_data("finalMAE_combo_matrix_small")
  res_l <- convert_combo_data_to_dt(test_mae[[1]])
  
  # expected assays converted
  exp_as <- as.character(get_combo_assay_names())
  expect_identical(sort(names(res_l)), sort(exp_as))
  
  # check content of data.table
  expect_true(nrow(res_l[[1]]) > 1 && ncol(res_l[[1]]) > 1)
  exp_idfs <- get_prettified_identifiers(c("drug_name", "drug_name2", "cellline"), simplify = FALSE)
  expect_true(all(exp_idfs %in% colnames(res_l[[1]])))
  
  # errors
  expect_error(
    convert_combo_data_to_dt(data.table::data.table(a = 1)),
    paste("Assertion on 'se' failed: Must inherit from class 'SummarizedExperiment',",
          "but has classes 'data.table','data.frame'."), 
    fixed = TRUE
  )
  err_msg <-
    "Assertion on 'dummy' failed. Must be element(s) of {'excess, scores, "
  err_msg2 <-
    "isobolograms'} set."
  expect_error(
    convert_combo_data_to_dt(test_mae[[1]], c_assays = c("dummy", "excess")),
    sprintf("%s%s", err_msg, err_msg2),
    fixed = TRUE
  )
  expect_error(
    convert_combo_data_to_dt(test_mae[[1]], normalization_type = 1),
    "Assertion on 'normalization_type' failed: Must be of type 'character', not 'double'.",
    fixed = TRUE
  )
  expect_error(
    convert_combo_data_to_dt(test_mae[[1]], prettify = "true"),
    "Assertion on 'prettify' failed: Must be of type 'logical flag', not 'character'.",
    fixed = TRUE
  )
})

test_that("shorten_normalization_type_name", {
  ### expected values
  expect_identical("GR", shorten_normalization_type_name("GRvalue"))
  
  ### errors
  err_msg <- "Assertion on 'x' failed: Must be element of set"
  expect_error(shorten_normalization_type_name("invalid"), err_msg)
})

test_that("define_matrix_grid_positions", {
  conc <- c(0, 10^(seq(-3, 1, 0.5)))
  conc_2 <- rep(conc[2:3], length.out = NROW(conc))

  output <- data.table::data.table(
    conc = round_concentration(conc),
    log10conc = log10(round_concentration(conc)),
    pos = log10(round_concentration(conc)),
    marks = sprintf("%.2g", conc)
  )
  output$pos[1] <- 
    2 * log10(round_concentration(conc))[2] - log10(round_concentration(conc))[3] - log10(1.5)
  
  res <- define_matrix_grid_positions(conc, conc_2)
  expect_is(res, "list")
  expect_length(res, 2)
  expect_equal(names(res), c("axis_1", "axis_2"))
  expect_equal(dim(res[[1]]), c(NROW(conc), 4))
  expect_equal(dim(res[[2]]), c(NROW(unique(conc_2)), 4))
  expect_equal(res[[1]]$conc_1, round_concentration(conc))
  expect_equal(res[[1]]$log10conc_1, output$log10conc)
  expect_equal(res[[1]]$pos_y, output$pos)
  expect_equal(res[[1]]$marks_y, output$marks)
  
  res_2 <- define_matrix_grid_positions(conc, NA)
  expect_is(res_2, "list")
  expect_length(res_2, 2)
  expect_equal(names(res_2), c("axis_1", "axis_2"))
  expect_equal(res_2$axis_1, res$axis_1)
  expect_equal(dim(res_2$axis_2), c(0, 4))
  
  res_3 <-  define_matrix_grid_positions(conc, c(1.2))
  expect_is(res_3, "list")
  expect_length(res_3, 2)
  expect_equal(res_3[["axis_2"]]$pos_x, log10(1.2))
  
  res_4 <-  define_matrix_grid_positions(conc, c(0, 1.2))
  expect_is(res_4, "list")
  expect_length(res_4, 2)
  expect_equal(res_4[["axis_2"]][conc_2 == 0, ]$pos_x, log10(1.2) - 0.5)
  
  expect_error(define_matrix_grid_positions(conc, LETTERS[1:5]))
  expect_error(define_matrix_grid_positions(NULL, conc))
})

test_that("round_concentration", {
  x <- c(0.00175, 0.00324, 0.0091)
  expect_equal(round_concentration(x), c(0.00175, 0.00324, 0.00910))
  expect_equal(round_concentration(x, ndigit = 1), c(0.002, 0.003, 0.010))
  
  expect_error(round_concentration(LETTERS[1:5]))
  expect_error(round_concentration(NULL))
  expect_error(round_concentration(x, ndigit = 1.5))
  expect_error(round_concentration(x, ndigit = "str"))
})
