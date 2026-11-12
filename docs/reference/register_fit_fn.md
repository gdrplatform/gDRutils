# Register or update a fit function

Adds a fit function to the registry read by
[`get_fit_fns`](https://gdrplatform.github.io/gDRstyle/reference/get_fit_fns.md),
or overwrites an entry of the same name. The entry is available
immediately in the current R session and is not saved anywhere: register
it again after restarting R, or from the registering package's
`.onLoad()`.

## Usage

``` r
register_fit_fn(name, package, fn_name = name, data_type, available = NULL)
```

## Arguments

- name:

  string; unique identifier, conventionally the function name.

- package:

  string; package exporting the function.

- fn_name:

  string; the exported name, when it differs from `name`.

- data_type:

  character vector of fit profile names the function is written for,
  e.g. `"single-agent"`. A vector rather than a string: a function may
  serve more than one profile.

- available:

  function of no arguments returning a flag, or `NULL`. Evaluated by
  [`get_fit_fns`](https://gdrplatform.github.io/gDRstyle/reference/get_fit_fns.md)
  to decide whether the function can run here — `musyc_available` for a
  fit needing Python, for instance. `NULL` means available whenever the
  package is installed.

## Value

Invisibly, the registered entry.

## Details

Registering a function does not change how it is called.
[`apply_fit`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.md)
takes the function itself, with or without an entry here.

## See also

[`get_fit_fns`](https://gdrplatform.github.io/gDRstyle/reference/get_fit_fns.md),
[`register_fit_profile`](https://gdrplatform.github.io/gDRstyle/reference/register_fit_profile.md)

## Examples

``` r
register_fit_fn(
  "fit_drug_response_metrics",
  package = "gDRutils",
  data_type = "single-agent"
)
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
