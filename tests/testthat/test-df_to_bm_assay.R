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
  expect_identical(
    as.matrix(result_default),
    as.matrix(result_precomputed)
  )
})

test_that("df_to_bm_assay validates precomputed_metadata type", {
  dt <- data.table::data.table(Gnumber = "G1", clid = "CL1", x = 1)
  expect_error(
    df_to_bm_assay(dt, precomputed_metadata = "not_a_list"),
    "list"
  )
})
