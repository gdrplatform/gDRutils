# Logistic fit

Fit a logistic curve to drug response data.

## Usage

``` r
logisticFit(
  concs,
  norm_values,
  std_norm_values = NA,
  x_0 = 1,
  priors = NULL,
  lower = NULL,
  range_conc = c(0.005, 5),
  force_fit = FALSE,
  pcutoff = 0.05,
  cap = 0.1,
  n_point_cutoff = 4,
  capping_fold = 5
)
```

## Arguments

- concs:

  concentrations that have not been transformed into log space.

- norm_values:

  normalized response values (Untreated = 1).

- std_norm_values:

  std of values.

- x_0:

  upper limit. Defaults to `1`. For co-treatments, this value should be
  set to `NA`.

- priors:

  numeric vector containing starting values for all. mean parameters in
  the model. Overrules any self starter function.

- lower:

  numeric vector of lower limits for all parameters in a 4-param model.

- range_conc:

  range of concentration for calculating AOC_range.

- force_fit:

  boolean indicating whether or not to force a parameter-based fit.

- pcutoff:

  numeric of pvalue significance threshold above or equal to which to
  use a constant fit.

- cap:

  numeric value capping `norm_values` to stay below (`x_0` + cap).

- n_point_cutoff:

  integer indicating number of unique concentrations required to fit
  curve.

- capping_fold:

  Integer value of the fold number to use for capping IC50/GR50. Default
  is `5`.

## Value

data.table with metrics and fit parameters.

## Details

Implementation of the genedata approach for curve fit:
https://screener.genedata.com/documentation/display/DOC21/Business-Rules-for-Dose-Response-Curve-Fitting,-Model-Selection,-and-Fit-Validity.html
\#nolint

The output parameter names correspond to the following definitions:

- x_mean:

  The mean of a given dose-response metric

- x_AOC_range:

  The range of the area over the curve

- x_AOC:

  The area over the GR curve or, respectively, under the relative cell
  count curve, averaged over the range of concentration values

- xc50:

  The concentration at which the effect reaches a value of 0.5 based on
  interpolation of the fitted curve

- x_max:

  The maximum effect of the drug

- ec50:

  The drug concentration at half-maximal effect

- x_inf:

  The asymptotic value of the sigmoidal fit to the dose-response data as
  concentration goes to infinity

- x_0:

  The asymptotic metric value corresponding to a concentration of 0 for
  the primary drug

- h:

  The hill coefficient of the fitted curve, which reflects how steep the
  dose-response curve is

- r2:

  The goodness of the fit

- x_sd_avg:

  The standard deviation of GR/IC

- fit_type:

  This will be given by one of the following:

  - "DRC4pHillFitModel" Successfully fit with a 4-parameter model

  - "DRC3pHillFitModelFixS0" Successfully fit with a 3-parameter model

  - "DRCConstantFitResult" Successfully fit with a constant fit

  - "DRCTooFewPointsToFit" Not enough points to run a fit

  - "DRCInvalidFitResult" Fit was attempted but failed

- maxlog10Concentration:

  The highest log10 concentration

- N_conc:

  Number of unique concentrations

## Examples

``` r
logisticFit(
c(0.001, 0.00316227766016838, 0.01, 0.0316227766016838),
c(0.9999964000144, 0.999964001439942, 0.999640143942423, 0.996414342629482),
rep(0.1, 4),
priors = c(2, 0.4, 1, 0.00658113883008419)
)
#>       x_mean        x_AOC x_AOC_range  xc50     x_max      ec50     x_inf   x_0
#>        <num>        <num>       <num> <num>     <num>     <num>     <num> <num>
#> 1: 0.9993262 0.0006738196 0.003303951   Inf 0.9964143 0.0187753 0.9958965     1
#>           h        r2      p_value          rss maxlog10Concentration N_conc
#>       <num>     <num>        <num>        <num>                 <num>  <int>
#> 1: 3.711073 0.9998954 0.0001046336 9.435296e-10                  -1.5      4
#>    x_sd_avg               fit_type
#>       <num>                 <char>
#> 1:      0.1 DRC3pHillFitModelFixS0
```
