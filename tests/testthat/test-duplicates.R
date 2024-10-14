test_that("has_dt_duplicated_rows works as expected", {
  
  dt_iris <- data.table::data.table(iris)
  expect_true(has_dt_duplicated_rows(dt_iris))
  expect_false(has_dt_duplicated_rows(dt_iris[1:100, ]))
  
  expect_true(has_dt_duplicated_rows(dt_iris[1:10, ], col_names = c("Sepal.Length", "Species")))
  expect_error(has_dt_duplicated_rows(iris), "Assertion on 'dt' failed")
  expect_error(
    has_dt_duplicated_rows(dt_iris, col_names = "invalid_value"),
    "Assertion on 'col_names' failed"
  )
  
})

test_that("get_duplicated_rows works as expected", {
  DF1co <- S4Vectors::DataFrame("Gnumber" = c("G0123456.1-1", "G0123456.2-2", "G1234567.1-1"),
                                "DrugName" = c("drug_name1", "drug_name1", "drug_name2"),
                                "Gnumber_2" = c("G9876543.1-1", "G9876543.1-1", "G9876543.1-1"),
                                "DrugName_2" = c("codrug_name1", "codrug_name1", "codrug_name1"),
                                "Concentration_2" = c("untreated", "untreated", "untreated"))
  
  
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
  # single column with non-default output
  expect_equal(
    get_duplicated_rows(DF1co, col_names = "DrugName", output = "data"),
    DF1co[1:2, ]
  )
  
  expect_error(get_duplicated_rows(DF1co, c("DrugName", "Fake Column")),
               "Assertion on 'all(col_names %in% colnames(x))' failed: Must be TRUE.", fixed = TRUE)
})

test_that("[has|get]_assay_dt_duplicated_rows works as expected", {
 
  # single-agent data
  sdata <- get_synthetic_data("finalMAE_small")
  smetrics_data <- convert_se_assay_to_dt(sdata[[get_supported_experiments("sa")]], "Metrics")
  smetrics_dup_data <- rbind(smetrics_data, smetrics_data[1:10, ])
  expect_false(has_assay_dt_duplicated_rows(smetrics_data))
  expect_true(has_assay_dt_duplicated_rows(smetrics_dup_data))
  
  expect_equal(get_assay_dt_duplicated_rows(smetrics_dup_data),
               c(1:10, 201:210))
  expect_equal(
    rbind(smetrics_data[1:10, ], smetrics_data[1:10, ]),
    get_assay_dt_duplicated_rows(smetrics_dup_data, output = "data")
  )
  expect_equal(get_assay_dt_duplicated_rows(smetrics_data), integer(0))
  empty_dt <- get_assay_dt_duplicated_rows(smetrics_data, output = "data")
  expect_true(nrow(empty_dt) == 0)
  expect_is(empty_dt, "data.table")
  
  # combo data
  cdata <- get_synthetic_data("finalMAE_combo_matrix_small")
  cscores_data <- convert_se_assay_to_dt(cdata[[1]], "scores")
  cscores_dup_data <- rbind(cscores_data, cscores_data[1:10, ])
  expect_false(has_assay_dt_duplicated_rows(cscores_data))
  expect_true(has_assay_dt_duplicated_rows(cscores_dup_data))
  
  expect_equal(get_assay_dt_duplicated_rows(cscores_dup_data),
               c(1:10, 25:34))
  expect_equal(
    rbind(cscores_data[1:10, ], cscores_data[1:10, ]),
    get_assay_dt_duplicated_rows(cscores_dup_data, output = "data")
  )
  expect_equal(get_assay_dt_duplicated_rows(cscores_data), integer(0))
  empty_dt <- get_assay_dt_duplicated_rows(cscores_data, output = "data")
  expect_true(nrow(empty_dt) == 0)
  expect_is(empty_dt, "data.table")
})

test_that("throw_msg_if_duplicates works as expected", {
 
  sdata <- get_synthetic_data("finalMAE_small")
  smetrics_data <- convert_se_assay_to_dt(sdata[[1]], "Metrics")
  smetrics_dup_data <- rbind(smetrics_data, smetrics_data[1:10, ])
  
  exp_msg <- "rows are duplicated"
  expect_error(throw_msg_if_duplicates(smetrics_dup_data, "Metrics"), exp_msg)
  expect_warning(throw_msg_if_duplicates(smetrics_dup_data, "Metrics", msg_f = warning), exp_msg)
})
