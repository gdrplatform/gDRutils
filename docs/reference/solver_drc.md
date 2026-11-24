# Hill fit solver using drc

The default solver of
[`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.md):
[`drc::drm()`](https://rdrr.io/pkg/drc/man/drm.html) with
[`drc::LL.3u()`](https://rdrr.io/pkg/drc/man/LL.3.html) (upper asymptote
fixed) or [`drc::LL.4()`](https://rdrr.io/pkg/drc/man/LL.4.html) (upper
asymptote fitted), within the bounds of the model. With bounds, drc
optimises with `stats::optim(method = "L-BFGS-B")`.

## Usage

``` r
solver_drc(conc, x, model, control = list())
```

## Arguments

- conc:

  numeric vector of concentrations (not logged).

- x:

  numeric vector of responses, the same length as `conc`.

- model:

  list describing the model to fit: `free` parameter names, `fixed`
  values (`x_0` for the 3-parameter model), and `start`, `lower`,
  `upper` over the free parameters.

- control:

  named list of [`drc::drmc()`](https://rdrr.io/pkg/drc/man/drmc.html)
  settings overriding the defaults (`relTol = 1e-4`, `errorm = FALSE`,
  `noMessage = TRUE`, `rmNA = TRUE`). With the bounds of the model drc
  optimises with L-BFGS-B, which ignores `relTol`.

## Value

list with `params` (named `h`, `x_inf`, `x_0`, `ec50`), `rss`,
`converged` and `message`.

## Details

A solver is any function with this signature returning the same list,
and can be passed to `fit_drug_response_metrics(solver = )` or through
`apply_fit(fit_args = list(solver = ))`.

## Examples

``` r
model <- list(free = c("h", "x_inf", "ec50"), fixed = list(x_0 = 1),
              start = c(h = 2, x_inf = 0.4, ec50 = 0.1),
              lower = c(h = 0.1, x_inf = 0, ec50 = 1e-4),
              upper = c(h = 5, x_inf = 1, ec50 = 100))
solver_drc(c(0.001, 0.01, 0.1, 1, 10), c(0.95, 0.8, 0.5, 0.2, 0.1), model)
#> $params
#>          h      x_inf        x_0       ec50 
#> 0.65230324 0.05710528 1.00000000 0.07913255 
#> 
#> $rss
#> [1] 0.00018015
#> 
#> $converged
#> [1] TRUE
#> 
#> $message
#> [1] NA
#> 
```
