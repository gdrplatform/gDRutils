#library(testthat); library(gDRutils)

test_that("has_nested_field works as expected", {
  m <- 10
  n <- 5
  df <- S4Vectors::DataFrame(clids = rep(LETTERS, each = m),
                   drugs = rep(letters[1:n], length(LETTERS) * m/n),
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
                   drugs = rep(letters[1:n], length(LETTERS) * m/n),
                   drug_name = paste0(rep(letters[1:n], length(LETTERS) * m/n), "A"),
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
  expect_equal(dim(obs), c(n, length(LETTERS)))
})

test_that("aggregate_assay works as expected", {
  m <- 10
  n <- 5
  df <- S4Vectors::DataFrame(clids = rep(LETTERS, each = m),
                   drugs = rep(letters[1:n], length(LETTERS) * m/n),
                   group = rep(LETTERS, each = m),
                   GR_50 = rep(seq(length(LETTERS)), each = m),
                   IC_50 = rep(seq(m), length(LETTERS))
  )
  asy <- BumpyMatrix::splitAsBumpyMatrix(df[, c("group", "GR_50", "IC_50")], row = df$drugs, column = df$clids)
   
  expect_error(aggregate_assay(asy, FUN = sum, by = c("clids")), regexp = "specified 'by' columns: 'clids' are not present in 'asy'")

  obs_asy <- aggregate_assay(asy, FUN = sum, by = c("group"))
  obs_df <- BumpyMatrix::unsplitAsDataFrame(obs_asy)

  expect_true(is(obs_asy, "BumpyMatrix"))
  expect_equal(rownames(obs_asy), rownames(asy))
  expect_equal(colnames(obs_asy), colnames(asy))
  expect_equal(nrow(obs_df), length(LETTERS) * n)
  expect_true(all(obs_df$column == obs_df$group))
})
