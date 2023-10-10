test_that("flatten works as expected", {
  n <- 4
  m <- 5
  grid <- expand.grid(normalization_type = c("GR", "RV"),
                      source = c("GDS", "GDR"))
  repgrid <- do.call("rbind", rep(list(grid), m))
  repgrid$wide <- seq(m * n)
  repgrid$id <- rep(LETTERS[1:m], each = n)
  repgrid$id2 <- rep(paste0(LETTERS[1:m], "2"), each = n)
  repgrid$constant <- "constant"
  repgrid <- data.table::setDT(repgrid)

  groups <- colnames(grid)
  wide_cols <- c("wide")

  # data.frame
  out <- flatten(repgrid, groups = groups, wide_cols = wide_cols)
  expect_equal(dim(out), c(m, n * length(wide_cols) + length(setdiff(colnames(repgrid), c(groups, wide_cols)))))
  expect_setequal(colnames(out),
                  c("id", "id2", "constant", "GR_GDS_wide",
                    "RV_GDS_wide", "GR_GDR_wide", "RV_GDR_wide"))

  # data.table
  repgrid2 <- data.table::as.data.table(repgrid)
  out2 <- flatten(repgrid2, groups = groups, wide_cols = wide_cols)
  expect_equal(dim(out2), c(m, n * length(wide_cols) + length(setdiff(colnames(repgrid2), c(groups, wide_cols)))))
  expect_setequal(colnames(out2),
                  c("id", "id2", "constant", "GR_GDS_wide",
                    "RV_GDS_wide", "GR_GDR_wide", "RV_GDR_wide"))

  # Remove one fit_source x normalization_type combination.
  repgrid3 <- repgrid[!(repgrid$normalization_type == "GR" & repgrid$source == "GDS"), ]
  out3 <- flatten(repgrid3, groups = groups, wide_cols = wide_cols)
  expect_equal(dim(out3), c(m, (n - 1) * length(wide_cols) + length(setdiff(colnames(repgrid3), c(groups, wide_cols)))))
  expect_setequal(colnames(out3),
                  c("id", "id2", "constant",
                    "RV_GDS_wide", "GR_GDR_wide", "RV_GDR_wide"))
  
  repgrid3$x <- runif(nrow(repgrid3))
  out4 <- flatten(repgrid3, groups = "normalization_type", wide_cols = "x")
  expect_setequal(colnames(out4),
                  c("source", "wide", "id", "id2",
                    "constant", "GRvalue", "RelativeViability"
                  ))
})
