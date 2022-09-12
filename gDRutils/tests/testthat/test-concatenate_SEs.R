#library(testthat); library(gDRutils)

test_that("aggregate_assay works as expected", {
  m <- 10
  n <- 5
  df <- S4Vectors::DataFrame(clids = rep(LETTERS, each = m),
                   drugs = rep(letters[1:n], length(LETTERS)*m/n),
                   group = rep(LETTERS, each = m),
                   GR_50 = rep(seq(length(LETTERS)), each = m),
                   IC_50 = rep(seq(m), length(LETTERS))
  )
  asy <- BumpyMatrix::splitAsBumpyMatrix(df[, c("group", "GR_50", "IC_50")], row = df$drugs, column = df$clids)
   
  expect_error(aggregate_assay(asy, aggregation_fxn = sum, by = c("clids")), regexp = "specified 'by' columns: 'clids' are not present in 'asy'")

  obs_asy <- aggregate_assay(asy, aggregation_fxn = sum, by = c("group"))
  obs_df <- BumpyMatrix::unsplitAsDataFrame(obs_asy)
  expect_true(is(obs_asy, "BumpyMatrix"))
  expect_equal(rownames(obs_asy), rownames(asy))
  expect_equal(colnames(obs_asy), colnames(asy))
  expect_equal(nrow(obs_df), length(LETTERS) * n)
  expect_true(all(obs_df$column == obs_df$group))
})
