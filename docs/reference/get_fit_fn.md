# Get a registered fit function

Resolves a registry entry to the function itself.

## Usage

``` r
get_fit_fn(name)
```

## Arguments

- name:

  string; the registered identifier.

## Value

The fit function. Stops if the name is not registered, or if its package
is not installed.

## See also

[`get_fit_fns`](https://gdrplatform.github.io/gDRstyle/reference/get_fit_fns.md)

## Examples

``` r
get_fit_fn("fit_drug_response_metrics")
#> function (avg_dt, x_col = "x", capping_fold = 5, range_conc = c(0.005, 
#>     5), pcutoff = 0.05, n_point_cutoff = 4L, force_fit = FALSE, 
#>     cap = 0.1) 
#> {
#>     .fit_drug_response_metrics_impl(avg_dt, x_0 = 1, x_col = x_col, 
#>         capping_fold = capping_fold, range_conc = range_conc, 
#>         pcutoff = pcutoff, n_point_cutoff = n_point_cutoff, force_fit = force_fit, 
#>         cap = cap)
#> }
#> <bytecode: 0x5690d69607c8>
#> <environment: namespace:gDRutils>
#> attr(,"fit_fn_info")
#> attr(,"fit_fn_info")$model
#> [1] "3-parameter log-logistic (drc::LL.3u), upper asymptote fixed at 1"
#> 
#> attr(,"fit_fn_info")$description
#> [1] "Hill curve of the response against concentration, fitted per normalization type. Falls back to a constant fit when fewer than n_point_cutoff concentrations are available or, unless force_fit = TRUE, when the F-test p-value is at least pcutoff."
#> 
#> attr(,"fit_fn_info")$reference
#> [1] "Ritz et al. 2015, PLOS ONE, doi:10.1371/journal.pone.0146021"
#> 
```
