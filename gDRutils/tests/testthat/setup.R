# Set up synthetic data for testing.

.create_data <- function(conc, params, params_GR) {
  data.frame(
    Concentration = conc,
    std_RelativeViability = 0.1,
    std_GRvalue = 0.1,
    RelativeViability =
      logistic_4parameters(
        conc, 
        params$x_inf, 
        params$x_0,
        params$ec50, 
        params$h
      ),
    GRvalue =
      logistic_4parameters(
        conc,
        params_GR$x_inf,
        params_GR$x_0,
        params_GR$ec50,
        params_GR$h
      )
  )
}

params <- data.frame(h = 2, x_inf = 0.1, x_0 = 1, ec50 = 0.5)
params_GR <- params
params_GR$x_inf <- -0.4

expected <- rbind(params, params_GR)
rownames(expected) <- c("RV_gDR", "GR_gDR")
expected_dims <- c(2, 16)

conc <- 10 ^ (seq(-3, 1, 0.5))

# Normal curve.
df_resp <- .create_data(conc, params, params_GR) 

# Above 0.5 curve.
params_above <- data.frame(h = 2, x_inf = 0.5, x_0 = 1, ec50 = 0.75)
df_resp_above <- .create_data(conc, params_above, params_above)

# Below 0.5 curve. 
params_below <- data.frame(h = 2, x_inf = 0, x_0 = 0.4, ec50 = 0.2)
df_resp_below <- .create_data(conc, params_below, params_below)

n <- 64
md_df <- data.frame(
  Gnumber = rep(c("vehicle", "untreated", paste0("G", seq(2))), each = 16), 
  DrugName = rep(c("vehicle", "untreated", paste0("GN", seq(2))), each = 16), 
  clid = paste0("C", rep_len(seq(4), n)),
  CellLineName = paste0("N", rep_len(seq(4), n)),
  replicates = rep_len(paste0("R", rep(seq(4), each = 4)), 64),
  drug_moa = "inhibitor",
  ReferenceDivisionTime = rep_len(c(120, 60), n),
  Tissue = "Lung",
  parental_identifier = "CL12345",
  Duration = 160
)

data_df <- data.frame(
  Concentration = rep(c(0, 0, 1, 3), each = 16), 
  ReadoutValue = runif(n, 1000, 5000),
  BackgroundValue = runif(n, 0, 1),
  WellRow = rep_len(LETTERS[1:8], n),
  WellColumn = rep_len(seq(3), n),
  experimenter = "Bob Ross"
)

test_df <- cbind(md_df, data_df)
