# The pipeline's fit, with the concentration at further response levels

Runs
[`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md)
and appends the concentration at each requested inhibition level —
`xc90` for IC90, `xc75` for IC75, and so on.

## Usage

``` r
fit_drug_response_metrics_xcx(
  avg_dt,
  x_col = "x",
  levels = c(50, 90),
  capping_fold = 5,
  range_conc = c(0.005, 5),
  pcutoff = 0.05,
  n_point_cutoff = 4L,
  force_fit = FALSE,
  cap = 0.1
)
```

## Arguments

- avg_dt:

  `data.table` for one drug, cell line and normalization type, as
  [`apply_fit`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.md)
  supplies it.

- x_col:

  string; the response column.

- levels:

  numeric vector of inhibition percentages between 1 and 99.

- capping_fold, range_conc, pcutoff, n_point_cutoff, force_fit, cap:

  passed to
  [`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md)
  unchanged.

## Value

The named list
[`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md)
returns, plus one element per requested level named `xc<level>`.

## Details

The curve is not changed. This calls the pipeline's own fit with the
arguments it is given, then inverts the resulting model with
[`predict_conc_from_efficacy`](https://gdrplatform.github.io/gDRstyle/reference/predict_conc_from_efficacy.md)
at the requested levels. Use it when the fit is right and only `xc50` is
too few points to read off it; use one of the alternative curve models
when the *shape* is wrong.

**`levels` are inhibition percentages.** `levels = 90` asks for the
concentration at which the response reaches `0.1`, not `0.9`. This
matches `CIScore_80` on the combination side. `xc50` is consistent with
either reading only because 50 is its own complement, so the distinction
matters at every other level.

`xc50` is never recomputed: it is already produced by the fit and means
the same thing, so asking for level 50 returns the fit's own value
rather than a second column of the same name.

The `xc` prefix is deliberate and stays normalization-agnostic, so the
same column name is honest on the GR scale, where "IC" is not.

Levels that the curve never reaches come back as `Inf`, and levels that
extrapolate more than `capping_fold` beyond the tested concentrations
are capped to `±Inf` by
[`cap_xc50`](https://gdrplatform.github.io/gDRstyle/reference/cap_xc50.md)
— the same treatment `xc50` already receives, not a parallel rule.

The fit-control arguments are passed through so a dataset processed with
non-default settings can be matched exactly. `fit_SE()` records those
settings in the provenance recipe, so they can be read back with
[`get_fit_provenance`](https://gdrplatform.github.io/gDRstyle/reference/fit_provenance.md)
and supplied here.

## Not registered for averaging

The new columns are concentrations, but they are not in
`get_header("metric_average_fields")$geometric_mean`. Do not average or
merge them as metrics:
[`average_biological_replicates_dt`](https://gdrplatform.github.io/gDRstyle/reference/average_biological_replicates_dt.md)
derives its field list from every numeric column and would take an
arithmetic mean of a quantity spanning orders of magnitude. Nothing
inside gDR averages a custom assay, and
[`merge_SE()`](https://gdrplatform.github.io/gDRstyle/reference/merge_SE.md)
fills with `NA` rather than averaging, so this is a limitation to
respect rather than a trap already sprung.

## See also

[`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md),
[`predict_conc_from_efficacy`](https://gdrplatform.github.io/gDRstyle/reference/predict_conc_from_efficacy.md),
[`cap_xc50`](https://gdrplatform.github.io/gDRstyle/reference/cap_xc50.md)

## Examples

``` r
dt <- data.table::data.table(
  Concentration = 10^(-3:3),
  x = c(1, .98, .9, .5, .2, .1, .08),
  normalization_type = "RV"
)
fit_drug_response_metrics_xcx(dt, levels = c(75, 90))
#> $normalization_type
#> [1] "RV"
#> 
#> $x_mean
#> [1] 0.5346938
#> 
#> $x_AOC
#> [1] 0.4653062
#> 
#> $x_AOC_range
#> [1] 0.2596304
#> 
#> $x_max
#> [1] 0.08
#> 
#> $x_sd_avg
#> [1] NaN
#> 
#> $N_conc
#> [1] 7
#> 
#> $maxlog10Concentration
#> [1] 3
#> 
#> $ec50
#> [1] 0.8603323
#> 
#> $xc50
#> [1] 1.070997
#> 
#> $h
#> h:(Intercept) 
#>     0.8944381 
#> 
#> $r2
#> [1] 0.9990855
#> 
#> $rss
#> [1] 0.0009667954
#> 
#> $p_value
#> [1] 2.529311e-08
#> 
#> $x_0
#> [1] 1
#> 
#> $x_inf
#> [1] 0.08895516
#> 
#> $fit_type
#> [1] "DRC3pHillFitModelFixS0"
#> 
#> $xc75
#> [1] 4.804318
#> 
#> $xc90
#> [1] 117.8409
#> 
```
