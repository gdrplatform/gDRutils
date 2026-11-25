# Check that recorded fits can be reproduced

Re-run each recorded fit with
[`refit_from_provenance`](https://gdrplatform.github.io/gDRstyle/reference/refit_from_provenance.md)
and compare the new rows with the stored ones.

## Usage

``` r
check_reproducibility(
  se,
  assay = NULL,
  fit_source = NULL,
  tolerance = sqrt(.Machine$double.eps),
  allow_source = FALSE,
  fit_id = NULL
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

- assay:

  string or `NULL`; check only fits written to this assay.

- fit_source:

  string or `NULL`; check only this fit source.

- tolerance:

  numeric; relative tolerance for numeric columns, as in
  [`all.equal`](https://rdrr.io/r/base/all.equal.html).

- allow_source:

  logical; evaluate fit functions stored as source code.

- fit_id:

  string or `NULL`; check only this recipe.

## Value

`data.table` with one row per checked fit: `assay`, `fit_source`,
`fit_id`, `reproducible` (`NA` when the fit could not be re-run),
`n_rows` and `detail`.

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[["single-agent"]]
se <- apply_fit(
  se, fit_drug_response_metrics, "single-agent",
  output_assay = "custom_hill", fit_source = "hill"
)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
check_reproducibility(se)
#>          assay fit_source       fit_id reproducible n_rows
#>         <char>     <char>       <char>       <lgcl>  <int>
#> 1:     Metrics        gDR          gDR        FALSE    200
#> 2: custom_hill       hill 74e28759241a         TRUE    200
#>                                                                                                                                                                                                                    detail
#>                                                                                                                                                                                                                    <char>
#> 1: Column 'ec50': Mean relative difference: 2.838961e-07; gDRutils::fit_drug_response_metrics was recorded with version 1.11.22 but 1.13.4 is installed; results may differ; the fit now has fit_id 7207b9e81092, not gDR
#> 2:                                                                                                                                                                                                          results match
```
