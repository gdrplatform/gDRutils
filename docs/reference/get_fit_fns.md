# List the registered fit functions

The fit functions known to this session, with the description each one
carries, for presenting a choice.

## Usage

``` r
get_fit_fns(data_type = NULL)
```

## Arguments

- data_type:

  character vector; keep only functions written for these profiles.
  `NULL` keeps all.

## Value

`data.frame`, one row per registered function, with columns `name`,
`package`, `fn_name`, `data_type` (comma separated), `model`,
`description`, `reference` and `available`. Zero rows when nothing
matches.

## Details

A function is listed as unavailable when its package is not installed,
or when its `available` predicate returns `FALSE` — a fit needing Python
is registered but not offered on a machine without it.

`model`, `description` and `reference` are read from the attributes
[`set_fit_fn_info`](https://gdrplatform.github.io/gDRstyle/reference/set_fit_fn_info.md)
attaches, and are `NA` for a function that does not carry them.

## Providers register when they load

An entry exists only once the package registering it has been loaded, so
a fit function from another package is missing from this table until
then. Call [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html)
on the packages you care about first if you need the full catalogue.

## See also

[`register_fit_fn`](https://gdrplatform.github.io/gDRstyle/reference/register_fit_fn.md),
[`get_fit_profiles`](https://gdrplatform.github.io/gDRstyle/reference/get_fit_profiles.md)

## Examples

``` r
get_fit_fns()
#>                           name  package                      fn_name
#> 1    fit_drug_response_metrics gDRutils    fit_drug_response_metrics
#> 2 fit_drug_response_metrics_4p gDRutils fit_drug_response_metrics_4p
#>      data_type
#> 1 single-agent
#> 2 single-agent
#>                                                               model
#> 1 3-parameter log-logistic (drc::LL.3u), upper asymptote fixed at 1
#> 2      4-parameter log-logistic (drc::LL.4), upper asymptote fitted
#>                                                                                                                                                                                                                                           description
#> 1 Hill curve of the response against concentration, fitted per normalization type. Falls back to a constant fit when fewer than n_point_cutoff concentrations are available or, unless force_fit = TRUE, when the F-test p-value is at least pcutoff.
#> 2                                                                         Hill curve of the response against concentration with a free upper asymptote, fitted per normalization type. Falls back to a constant fit like fit_drug_response_metrics().
#>                                                      reference available
#> 1 Ritz et al. 2015, PLOS ONE, doi:10.1371/journal.pone.0146021      TRUE
#> 2 Ritz et al. 2015, PLOS ONE, doi:10.1371/journal.pone.0146021      TRUE
get_fit_fns(data_type = "single-agent")
#>                           name  package                      fn_name
#> 1    fit_drug_response_metrics gDRutils    fit_drug_response_metrics
#> 2 fit_drug_response_metrics_4p gDRutils fit_drug_response_metrics_4p
#>      data_type
#> 1 single-agent
#> 2 single-agent
#>                                                               model
#> 1 3-parameter log-logistic (drc::LL.3u), upper asymptote fixed at 1
#> 2      4-parameter log-logistic (drc::LL.4), upper asymptote fitted
#>                                                                                                                                                                                                                                           description
#> 1 Hill curve of the response against concentration, fitted per normalization type. Falls back to a constant fit when fewer than n_point_cutoff concentrations are available or, unless force_fit = TRUE, when the F-test p-value is at least pcutoff.
#> 2                                                                         Hill curve of the response against concentration with a free upper asymptote, fitted per normalization type. Falls back to a constant fit like fit_drug_response_metrics().
#>                                                      reference available
#> 1 Ritz et al. 2015, PLOS ONE, doi:10.1371/journal.pone.0146021      TRUE
#> 2 Ritz et al. 2015, PLOS ONE, doi:10.1371/journal.pone.0146021      TRUE
```
