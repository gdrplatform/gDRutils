# Defunct: fit curves with the old fitting path

`fit_curves()` and `logisticFit()` were a second copy of the gDR Hill
fit with hard-coded priors. Every fit now goes through
[`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md),
whose priors and bounds come from `inst/extdata/fit_profiles.json` and
whose optimiser is a parameter
([`solver_drc`](https://gdrplatform.github.io/gDRstyle/reference/solver_drc.md)
by default).

## Usage

``` r
fit_curves(...)

logisticFit(...)
```

## Arguments

- ...:

  ignored.

## Value

Does not return; stops with a message pointing to
[`fit_drug_response_metrics()`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md).
