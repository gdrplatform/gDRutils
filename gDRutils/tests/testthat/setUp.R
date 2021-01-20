# Set up synthetic data for testing.

params <- data.frame(h = 2, x_inf = 0.1, x_0 = 1, c50 = 0.5)

params_GR <- params
params_GR$x_inf <- -0.4

conc <- 10 ^ (seq(-3, 1, 0.5))

df_resp <- data.frame(Concentration = conc, 
		      std_RelativeViability = 0.1,
		      std_GRvalue = 0.1,
		      RelativeViability =
			logistic_4parameters(conc, params$x_inf, params$x_0, 
			  params$c50, params$h),
		      GRvalue = 
                        logistic_4parameters(conc, params_GR$x_inf, 
		          params_GR$x_0, params_GR$c50, params_GR$h)
)

