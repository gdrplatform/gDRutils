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
        params$c50, 
        params$h
      ),
    GRvalue =
      logistic_4parameters(
        conc,
        params_GR$x_inf,
        params_GR$x_0,
        params_GR$c50,
        params_GR$h
      )
  )
}

params <- data.frame(h = 2, x_inf = 0.1, x_0 = 1, c50 = 0.5)
params_GR <- params
params_GR$x_inf <- -0.4

expected <- rbind(params, params_GR)
rownames(expected) <- c("RV_gDR", "GR_gDR")
expected_dims <- c(2, 16)

conc <- 10 ^ (seq(-3, 1, 0.5))

# Normal curve.
df_resp <- .create_data(conc, params, params_GR) 

# Above 0.5 curve.
params_above <- data.frame(h = 2, x_inf = 0.5, x_0 = 1, c50 = 0.75)
df_resp_above <- .create_data(conc, params_above, params_above)

# Below 0.5 curve. 
params_below <- data.frame(h = 2, x_inf = 0, x_0 = 0.4, c50 = 0.2)
df_resp_below <- .create_data(conc, params_below, params_below)
