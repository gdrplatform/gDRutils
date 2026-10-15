# Re-run a recorded fit

Rebuild the fit function from its recipe in `metadata(se)$provenance`
and fit the data again. The rows of `fit_source` in the output assay
(and in the summary assay, if any) are replaced; other rows are kept.

## Usage

``` r
refit_from_provenance(se, assay, fit_source, allow_source = FALSE)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

- assay:

  string; name of the assay holding the fit results.

- fit_source:

  string; value of the `fit_source` column.

- allow_source:

  logical; evaluate fit functions stored as source code.

## Value

`se` with the fit re-run.

## Details

A function from an installed package is taken from that package. A
warning is raised when the installed version or the code differs from
the recorded one, since results may then differ.

A function stored as source code is evaluated only with
`allow_source = TRUE`: the object may come from someone else, and
evaluating its code runs that code in your session.

## See also

[`check_reproducibility`](https://gdrplatform.github.io/gDRstyle/reference/check_reproducibility.md)

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
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
se <- refit_from_provenance(se, "custom_hill", "hill")
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
```
