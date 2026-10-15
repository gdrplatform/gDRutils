# Describe a fit function

Attach a short human-readable description to a fit function. It is
copied into the recipe that
[`apply_fit`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.md)
and
[`apply_fits`](https://gdrplatform.github.io/gDRstyle/reference/apply_fits.md)
record in `metadata(se)$provenance`, so a stored object says not only
which function produced a fit but also what that function does.

## Usage

``` r
set_fit_fn_info(
  fn,
  model = NULL,
  description = NULL,
  reference = NULL,
  env = NULL
)
```

## Arguments

- fn:

  the fit function to describe.

- model:

  string; name of the model, e.g. `"3-parameter log-logistic"`.

- description:

  string; what the function computes, in a sentence or two.

- reference:

  string; citation of the method, ideally with a DOI.

- env:

  optional function without arguments returning a named list of strings
  that describe the external runtime the fit depends on, e.g. the
  version of a Python package called through reticulate. It is called
  each time a fit is recorded.

## Value

`fn` with the description attached.

## See also

[`get_fit_provenance`](https://gdrplatform.github.io/gDRstyle/reference/fit_provenance.md)

## Examples

``` r
mean_fn <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
mean_fn <- set_fit_fn_info(
  mean_fn,
  model = "mean response",
  description = "Mean of the averaged response over all concentrations."
)
```
