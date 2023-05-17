test_that("has_nested_field works as expected", {
  m <- 10
  n <- 5
  df <- S4Vectors::DataFrame(clids = rep(LETTERS, each = m),
                             drugs = rep(letters[1:n], length(LETTERS) * m / n),
                             group = rep(LETTERS, each = m),
                             GR_50 = rep(seq(length(LETTERS)), each = m),
                             IC_50 = rep(seq(m), length(LETTERS))
  )
  asy <- BumpyMatrix::splitAsBumpyMatrix(df[, c("group", "GR_50", "IC_50")], row = df$drugs, column = df$clids)
  expect_true(gDRutils:::has_nested_field(asy, "GR_50"))
  expect_true(gDRutils:::has_nested_field(asy, c("GR_50", "group")))
  expect_false(gDRutils:::has_nested_field(asy, c("clids", "GR_50", "group")))
})


test_that(".transform_df_to_matrix works as expected", {
  m <- 10
  n <- 5
  df <- S4Vectors::DataFrame(clids = rep(LETTERS, each = m),
                             cellline_name = paste0(rep(LETTERS, each = m), "A"),
                             drugs = rep(letters[1:n], length(LETTERS) * m / n),
                             drug_name = paste0(rep(letters[1:n], length(LETTERS) * m / n), "A"),
                             group = rep(LETTERS, each = m),
                             GR_50 = rep(seq(length(LETTERS)), each = m),
                             IC_50 = rep(seq(m), length(LETTERS))
  )
  column_fields <- c("clids", "cellline_name")
  row_fields <- c("drugs", "drug_name")
  obs <- gDRutils:::.transform_df_to_matrix(df,
                                            row_fields = row_fields,
                                            column_fields = column_fields,
                                            nested_fields = c("group", "GR_50", "IC_50")
  )
  expect_equal(dim(obs$mat), c(n, length(LETTERS)))

  expect_error(gDRutils:::.transform_df_to_matrix(df, # TODO: FIX ME.
                                                  row_fields = row_fields,
                                                  column_fields = c(column_fields, row_fields[1]),
                                                  nested_fields = c("group", "GR_50", "IC_50")
  ))

})

test_that("demote_fields works as expected", {
  m <- 10
  n <- 5
  df <- S4Vectors::DataFrame(clids = rep(LETTERS, each = m),
                             cellline_name = paste0(rep(LETTERS, each = m), "A"),
                             drugs = rep(letters[1:n], length(LETTERS) * m / n),
                             drug_name = paste0(rep(letters[1:n], length(LETTERS) * m / n), "A"),
                             group = rep(LETTERS, each = m),
                             GR_50 = rep(seq(length(LETTERS)), each = m),
                             IC_50 = rep(seq(m), length(LETTERS))
  )

  # Demoting fields in rowData.
  column_fields <- c("clids", "cellline_name")
  row_fields <- c("drugs", "drug_name", "group")
  nested_fields <- c("GR_50", "IC_50")
  out <- gDRutils:::.transform_df_to_matrix(df,
                                            row_fields = row_fields,
                                            column_fields = column_fields,
                                            nested_fields = nested_fields
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list("test" = out$mat),
    rowData = out$rowData,
    colData = out$colData)

  obs <- demote_fields(se, "group")
  expect_equal(ncol(SummarizedExperiment::rowData(obs)), ncol(SummarizedExperiment::rowData(se)) - 1)
  expect_equal(colnames(SummarizedExperiment::assays(obs)[["test"]][1, 1][[1]]), c(nested_fields, "group"))
  expect_equal(nrow(obs), nrow(unique(df[, setdiff(row_fields, "group")])))
  expect_equal(ncol(obs), ncol(se))

  # Demoting fields in colData.
  column_fields <- c("clids", "cellline_name", "group")
  row_fields <- c("drugs", "drug_name")
  out <- gDRutils:::.transform_df_to_matrix(df,
                                            row_fields = row_fields,
                                            column_fields = column_fields,
                                            nested_fields = c("GR_50", "IC_50")
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list("test" = out$mat),
    rowData = out$rowData,
    colData = out$colData)
  obs <- demote_fields(se, "group")
  expect_equal(ncol(SummarizedExperiment::colData(obs)), ncol(SummarizedExperiment::colData(se)) - 1)
  expect_equal(colnames(SummarizedExperiment::assays(obs)[["test"]][1, 1][[1]]), c(nested_fields, "group"))
  expect_equal(ncol(obs), nrow(unique(df[, setdiff(column_fields, "group")])))
  expect_equal(nrow(obs), nrow(se))
})

test_that("promote_fields works as expected", {
  m <- 10
  n <- 5
  df <- S4Vectors::DataFrame(clids = rep(LETTERS, each = m),
                             cellline_name = paste0(rep(LETTERS, each = m), "A"),
                             drugs = rep(letters[1:n], length(LETTERS) * m / n),
                             drug_name = paste0(rep(letters[1:n], length(LETTERS) * m / n), "A"),
                             group = rep(LETTERS, each = m),
                             GR_50 = rep(seq(length(LETTERS)), each = m),
                             IC_50 = rep(seq(m), length(LETTERS))
  )
  column_fields <- c("clids", "cellline_name")
  row_fields <- c("drugs", "drug_name")
  nested_fields <- c("group", "GR_50", "IC_50")
  out <- gDRutils:::.transform_df_to_matrix(df,
                                            row_fields = row_fields,
                                            column_fields = column_fields,
                                            nested_fields = nested_fields
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list("test" = out$mat),
    rowData = out$rowData,
    colData = out$colData)

  obs <- promote_fields(se, "group", 1)
  expect_equal(nrow(obs), nrow(unique(df[, c("group", row_fields)])))
  expect_equal(colnames(SummarizedExperiment::rowData(obs)), c(row_fields, "group"))
  expect_equal(colnames(SummarizedExperiment::assays(obs)[["test"]][1, 1][[1]]), setdiff(nested_fields, "group"))
  expect_equal(ncol(obs), ncol(se))

  obs <- promote_fields(se, "group", 2)
  expect_equal(colnames(SummarizedExperiment::colData(obs)), c(column_fields, "group"))
  expect_equal(ncol(obs), nrow(unique(df[, c("group", column_fields)])))
  expect_equal(colnames(SummarizedExperiment::assays(obs)[["test"]][1, 1][[1]]), setdiff(nested_fields, "group"))
  expect_equal(nrow(obs), nrow(se))

  expect_error(promote_fields(se, "cellline_name", 2))
})

test_that("promote_fields and demote_fields are reversible operations", {
  m <- 10
  n <- 5
  df <- S4Vectors::DataFrame(clids = rep(LETTERS, each = m),
                             cellline_name = paste0(rep(LETTERS, each = m), "A"),
                             drugs = rep(letters[1:n], length(LETTERS) * m / n),
                             drug_name = paste0(rep(letters[1:n], length(LETTERS) * m / n), "A"),
                             group = rep(LETTERS, each = m),
                             GR_50 = rep(seq(length(LETTERS)), each = m),
                             IC_50 = rep(seq(m), length(LETTERS))
  )
  column_fields <- c("clids", "cellline_name")
  row_fields <- c("drugs", "drug_name")
  nested_fields <- c("group", "GR_50", "IC_50")
  out <- gDRutils:::.transform_df_to_matrix(df,
                                            row_fields = row_fields,
                                            column_fields = column_fields,
                                            nested_fields = nested_fields
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list("test" = out$mat),
    rowData = out$rowData,
    colData = out$colData)
  promoted_se <- promote_fields(se, "group", 1)
  demoted_promoted_se <- demote_fields(promoted_se, "group")
  obs <- S4Vectors::DataFrame(
    convert_se_assay_to_dt(demoted_promoted_se, "test", include_metadata = TRUE)
  )
  obs_df <- sort(obs[, setdiff(colnames(obs), c("rId", "cId"))])
  expect_equal(ncol(obs_df), ncol(df))
  expect_true(all(colnames(obs_df) %in% colnames(df)))
  expect_equal(sort(obs_df[, colnames(df)]), sort(df))
})

test_that("aggregate_assay works as expected", {
  m <- 10
  n <- 5
  df <- S4Vectors::DataFrame(clids = rep(LETTERS, each = m),
                             drugs = rep(letters[1:n], length(LETTERS) * m / n),
                             group = rep(LETTERS, each = m),
                             GR_50 = rep(seq(length(LETTERS)), each = m),
                             IC_50 = rep(seq(m), length(LETTERS))
  )
  asy <- BumpyMatrix::splitAsBumpyMatrix(df[, c("group", "GR_50", "IC_50")], row = df$drugs, column = df$clids)

  expect_error(aggregate_assay(asy, FUN = sum, by = c("clids")),
               regexp = "specified 'by' columns: 'clids' are not present in 'asy'")

  obs_asy <- aggregate_assay(asy, FUN = sum, by = c("group"))
  obs_df <- BumpyMatrix::unsplitAsDataFrame(obs_asy)

  expect_true(is(obs_asy, "BumpyMatrix"))
  expect_equal(rownames(obs_asy), rownames(asy))
  expect_equal(colnames(obs_asy), colnames(asy))
  expect_equal(nrow(obs_df), length(LETTERS) * n)
  expect_true(all(obs_df$column == obs_df$group))
})
