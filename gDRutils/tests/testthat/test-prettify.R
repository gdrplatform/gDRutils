test_that("prettify_flat_metrics works as expected", {
  x <- c("CellLineName", "Tissue",
         "Primary Tissue",
         "GR_gDR_x_mean", "GR_gDR_xc50", 
         "RV_GDS_x_mean", 
         "Concentration_2", "Gnumber_2", "Drug_3",
         "E_0", "GR_gDR_x_AOC_range"
  )
  
  y <- c("Gnumber", "Gnumber_2", "MyDrug", "MyDrug_2")

  obs <- prettify_flat_metrics(x, human_readable = FALSE)
  exp <- c("CellLineName", "Tissue",
           "Primary Tissue",
           "GR_mean", "GR50", 
           "GDS_RV_mean", 
           "Concentration_2", "Gnumber_2", "Drug_3",
           "E_0", "GR_AOC_range")
  expect_equal(obs, exp)

  # Human readable names work.
  obs <- prettify_flat_metrics(x, human_readable = TRUE)
  exp <- c("Cell Line Name", "Tissue",
           "Primary Tissue",
           "GR Mean", "GR50", 
           "RV Mean (GDS)",
           "Concentration 2", "Gnumber 2", "Drug 3",
           "E0", "GR AOC within set range")
  expect_equal(obs, exp)
  obs <- prettify_flat_metrics(y, human_readable = TRUE)
  exp <- c("Gnumber", "Gnumber 2", "My Drug", "My Drug 2")
  expect_equal(obs, exp)
})
