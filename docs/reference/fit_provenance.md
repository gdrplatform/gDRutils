# Get the provenance of fits stored in a SummarizedExperiment

[`apply_fit`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.md)
and
[`apply_fits`](https://gdrplatform.github.io/gDRstyle/reference/apply_fits.md)
record a recipe for every fit in `metadata(se)$provenance`. These
functions read it back.

## Usage

``` r
get_SE_provenance(se)

get_fit_provenance(se, assay, fit_source, data_source = NULL)

summarize_fit_provenance(se)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

- assay:

  string; name of the assay holding the fit results.

- fit_source:

  string; value of the `fit_source` column.

- data_source:

  string; for an SE built by
  [`merge_SE`](https://gdrplatform.github.io/gDRstyle/reference/merge_SE.md),
  the name of the merged object the rows come from (the `data_source`
  column). `NULL` reads the fits done on `se` itself.

## Value

`get_SE_provenance()`: the whole provenance entry, or `NULL` for objects
without one. `get_fit_provenance()`: one recipe.
`summarize_fit_provenance()`: a `data.table` with one row per recorded
fit.

## Details

A recipe contains:

- `fn`, `summary_fn`:

  record of the function: `name`, `package`, `version`, `git_sha` and
  `hash` for a function from an installed package; `source`, `captured`
  (values the function uses from its environment), `helpers` (source of
  helper functions it calls) and `dependencies` (versions of the
  packages it calls) for any other function; `info` set with
  [`set_fit_fn_info`](https://gdrplatform.github.io/gDRstyle/reference/set_fit_fn_info.md);
  `env` for declared runtime dependencies.

- `fit_args`:

  arguments passed to the fit function.

- `called_by`, `data_type`, `input_assay`, `output_assay`,
  `summary_assay`, `slicing_cols`, `slicing_values`, `merge`:

  settings of the call.

- `run`:

  R and gDRutils versions, hash of the input assay and the number of
  rows written.

Names listed under `not_stored` or `args_not_stored` were too large or
not plain data, so the fit cannot be fully reproduced from the object
alone.

For a summary assay, `get_fit_provenance()` returns the recipe of the
fit that wrote it. For an SE built by
[`merge_SE`](https://gdrplatform.github.io/gDRstyle/reference/merge_SE.md),
the recipes of the merged objects are kept per `data_source`.

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[["single-agent"]]
se <- apply_fit(
  se, fit_drug_response_metrics, "single-agent",
  output_assay = "custom_hill", fit_source = "hill_loose",
  fit_args = list(pcutoff = 0.1)
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
summarize_fit_provenance(se)
#>          assay fit_source                                            fn
#>         <char>     <char>                                        <char>
#> 1: custom_hill hill_loose gDRutils::fit_drug_response_metrics (1.11.19)
#>         fit_args
#>           <char>
#> 1: pcutoff = 0.1
#>                                                                model
#>                                                               <char>
#> 1: 3-parameter log-logistic (drc::LL.3u), upper asymptote fixed at 1
get_fit_provenance(se, "custom_hill", "hill_loose")$fit_args
#> $pcutoff
#> [1] 0.1
#> 
```
