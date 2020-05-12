library(testthat)

context("RVGRfits")

test_that("RVGRfits fails with expected errors", {
  expect_error(RVGRfits(list()),
  regexp = "inherits\\(df_, \"data.frame\"\\) is not TRUE")
})

test_that("Synthetic data", {

  # test 1: normal case
  params = c(2, .1, 1, .5)
  names(params) = c("h", "x_inf", "x_0", "c50")
  params_GR = params
  params_GR[2] = -.4
  conc = 10**(seq(-3,1,.5))

  df_resp = data.frame(Concentration = conc, std_RelativeViability = .1, std_GRvalue = .1,
    RelativeViability =
        logistic_4parameters(conc, params['x_inf'], params['x_0'], params['c50'], params['h']),
   GRvalue =
        logistic_4parameters(conc, params_GR['x_inf'], params_GR['x_0'],
          params_GR['c50'], params_GR['h'])
  )
  df_result = RVGRfits(df_resp)

  # tests: full fit
  apply(round(as.data.frame(df_result)[,names(params)],4) == rbind(params, params_GR),1,all)
  as.data.frame(df_result)[,'fit_type',drop=F] == 'DRC3pHillFitModelFixS0'

  df_result = RVGRfits(df_resp, e_0 = NA, GR_0 = NA)
  # tests: without the x_0 value
  apply(round(as.data.frame(df_result)[,names(params)],4) == rbind(params, params_GR),1,all)
  as.data.frame(df_result)[,'fit_type',drop=F] == 'DRC4pHillFitModel'

  df_resp
  df_result = RVGRfits(df_resp[3:6,])
  # tests:
  apply(is.na(as.data.frame(df_result[,names(params)])),1,all)
  as.data.frame(df_result)[,'fit_type',drop=F] == 'DRCTooFewPointsToFit'

}
