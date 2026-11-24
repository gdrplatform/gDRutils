# Deprecated: fit curves with the old interface

`fit_curves()` and `logisticFit()` were a second copy of the gDR Hill
fit, with hard-coded priors. They are now thin wrappers over
[`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md),
the single gDR Hill fit, whose priors and bounds come from
`inst/extdata/fit_profiles.json` and whose optimiser is a parameter
([`solver_drc`](https://gdrplatform.github.io/gDRstyle/reference/solver_drc.md)
by default). They are deprecated and will be defunct in the next
release: call `fit_drug_response_metrics(avg_dt, x_0 = )` once per
normalization type instead.

## Usage

``` r
fit_curves(
  df_,
  series_identifiers,
  e_0 = 1,
  GR_0 = 1,
  n_point_cutoff = 4,
  range_conc = c(0.005, 5),
  force_fit = FALSE,
  pcutoff = 0.05,
  cap = 0.1,
  normalization_type = c("GR", "RV")
)

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

- df\_:

  `data.table` or `DFrame` with the concentrations in
  `series_identifiers`, the responses in `x`, optionally their standard
  deviation in `x_std`, and `normalization_type`.

- series_identifiers:

  string, the concentration column of `df_`.

- e_0, GR_0:

  upper asymptote of the RV and of the GR fit.

- n_point_cutoff, range_conc, force_fit, pcutoff, cap, capping_fold:

  passed to
  [`fit_drug_response_metrics()`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md).

- normalization_type:

  normalization types to fit, among `"GR"` and `"RV"`.

- concs:

  concentrations, not logged.

- norm_values:

  responses, the same length as `concs`.

- std_norm_values:

  their standard deviation.

- x_0:

  upper asymptote, or `NA` to fit it.

- priors, lower:

  ignored: the start values and bounds are those of the normalization
  type's fit profile. A `lower` value below zero for `x_inf` still
  selects the GR profile, as the GR priors did.

## Value

`fit_curves()`: a `data.table` with one row per normalization type,
`fit_source` `"gDR"`, and rownames `"<type>_gDR"`. `logisticFit()`: a
one-row `data.table` of metrics and fit parameters.

## Examples

``` r
suppressWarnings(logisticFit(
  c(0.001, 0.00316227766016838, 0.01, 0.0316227766016838),
  c(0.9999964000144, 0.999964001439942, 0.999640143942423, 0.996414342629482),
  rep(0.1, 4)
))
#>    normalization_type    x_mean        x_AOC x_AOC_range     x_max x_sd_avg
#>                <char>     <num>        <num>       <num>     <num>    <num>
#> 1:                 RV 0.9994663 0.0005336782  0.07214082 0.9964143      0.1
#>    N_conc maxlog10Concentration      ec50  xc50        h    r2          rss
#>     <int>                 <num>     <num> <num>    <num> <num>        <num>
#> 1:      4                  -1.5 0.2030077   Inf 2.014503     1 3.425795e-13
#>         p_value   x_0     x_inf               fit_type
#>           <num> <num>     <num>                 <char>
#> 1: 3.799066e-08     1 0.8446057 DRC3pHillFitModelFixS0
```
