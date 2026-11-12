# Alternative dose-response curve models

Fit functions for
[`apply_fit`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.md),
each fitting one curve shape from `drc` to the averaged response of a
single-agent experiment.

## Usage

``` r
ll3_fit_fn(avg_dt, x_col = "x", conc_col = "Concentration")

ll5_fit_fn(avg_dt, x_col = "x", conc_col = "Concentration")

weibull1_fit_fn(avg_dt, x_col = "x", conc_col = "Concentration")

weibull2_fit_fn(avg_dt, x_col = "x", conc_col = "Concentration")

gompertz_fit_fn(avg_dt, x_col = "x", conc_col = "Concentration")

hormesis_fit_fn(avg_dt, x_col = "x", conc_col = "Concentration")
```

## Arguments

- avg_dt:

  `data.table` for one drug, cell line and normalization type, as
  [`apply_fit`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.md)
  supplies it.

- x_col:

  string; the response column.

- conc_col:

  string; the concentration column.

## Value

Named list of one row: the model's parameters, then `r2`, `converged`,
`ec50_extrapolated`, `n_conc` and `fit_type`.

## Details

These complement
[`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md),
which fits the 4-parameter log-logistic the processing pipeline uses.
Reach for one of these when that shape is the wrong one: `ll3_fit_fn`
when a curve is too sparse or too shallow to support four free
parameters, `ll5_fit_fn`, `weibull1_fit_fn`, `weibull2_fit_fn` and
`gompertz_fit_fn` when the response is asymmetric, and `hormesis_fit_fn`
for a low-dose stimulation that a log-logistic cannot represent at all.

Each returns the model's own parameters rather than the `Metrics`
schema. Parameters playing the same role across models share the names
[`logisticFit()`](https://gdrplatform.github.io/gDRstyle/reference/logisticFit.md)
uses — `h` for the slope, `x_0` and `x_inf` for the asymptotes, `ec50`
for the inflection point — and the rest are named for what they are.
Write them to an assay of their own: an assay of these is not the
`Metrics` schema, so nothing reading `Metrics` will understand it.

All of them also return `r2`, `converged`, `n_conc` and `fit_type`, so
two models fitted to the same data can be compared, and
`ec50_extrapolated`, which is `TRUE` when the fitted `ec50` lies outside
the concentrations tested. That happens whenever a drug is too weak to
reach its own midpoint on the plate: the fit can still be a good one,
but the `ec50` is not a measurement and should not be read as a potency.

A fit that does not converge returns a row of `NA` parameters with
`converged = FALSE`, rather than raising. A fit that failed is a result;
dropping the cell would make it indistinguishable from one that was
never fitted. `converged = FALSE` also covers a fit that finished on
parameters that say nothing about the data: one whose own model cannot
evaluate them, and one whose `ec50` ran a thousandfold outside the
tested concentrations, which happens on a response too flat to have a
midpoint.

`ec50` is bounded below at a positive value, because every model here
evaluates `log(dose/ec50)` and is undefined without it.

`ll3_fit_fn` and `hormesis_fit_fn` fix the lower asymptote at zero, so
neither suits a GR-normalized response that goes below it: the model
cannot reach the data and the misfit shows up as a poor `r2` rather than
as a failure.

## See also

[`apply_fit`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.md),
[`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md)

## Examples

``` r
dt <- data.frame(Concentration = 10^(-3:3), x = c(1, .98, .9, .5, .2, .1, .08))
ll3_fit_fn(dt)
#> $h
#> [1] 0.6366549
#> 
#> $x_0
#> [1] 1.032444
#> 
#> $ec50
#> [1] 1.127042
#> 
#> $r2
#> [1] 0.9900098
#> 
#> $converged
#> [1] TRUE
#> 
#> $ec50_extrapolated
#> [1] FALSE
#> 
#> $n_conc
#> [1] 7
#> 
#> $fit_type
#> [1] "DRC3pLogLogisticLowerFixed"
#> 
gompertz_fit_fn(dt)
#> $h
#> [1] NA
#> 
#> $x_inf
#> [1] NA
#> 
#> $x_0
#> [1] NA
#> 
#> $ec50
#> [1] NA
#> 
#> $r2
#> [1] NA
#> 
#> $converged
#> [1] FALSE
#> 
#> $ec50_extrapolated
#> [1] NA
#> 
#> $n_conc
#> [1] 7
#> 
#> $fit_type
#> [1] "DRC4pGompertz"
#> 
```
