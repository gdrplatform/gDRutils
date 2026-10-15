# Record provenance from code that does not use apply_fit()

[`apply_fit`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.md)
and
[`apply_fits`](https://gdrplatform.github.io/gDRstyle/reference/apply_fits.md)
record their recipes themselves. Code that writes results another way
uses these functions to add to the same `metadata(se)$provenance` entry.

## Usage

``` r
record_fit_provenance(
  se,
  assay,
  fit_source,
  fn,
  fit_args = list(),
  input_assay = NULL,
  replace = FALSE
)

record_processing_step(
  se,
  step,
  settings = list(),
  fns = list(),
  max_bytes = .provenance_max_value_bytes
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

- assay:

  string; assay the function wrote.

- fit_source:

  string; `fit_source` of the written rows. Assays without a
  `fit_source` column use it only as a key.

- fn:

  function of the whole SE that wrote `assay`, called as `fn(se, ...)`.

- fit_args:

  named list of the arguments `fn` was called with, other than `se`.

- input_assay:

  string or `NULL`; assay the results were computed from, hashed to tell
  later whether a re-run used the same data.

- replace:

  logical; `TRUE` when `fn` overwrote the whole assay, so the recipes of
  other fit sources in it no longer apply.

- step:

  string; name of the processing step.

- settings:

  named list of plain-data settings of the step.

- fns:

  named list of the functions the step used, e.g.
  `list(rate_fn = rate_fn)`.

- max_bytes:

  numeric; largest size in bytes of one setting that is stored. Larger
  settings are listed in `settings_not_stored`. Raise it for settings
  that are large by nature, such as the list of input files.

## Value

`se` with the provenance entry updated.

## Details

`record_fit_provenance()` is for a function of the whole SE that writes
an assay, e.g. a combination scoring step. The recipe keeps the function
and its arguments, so
[`refit_from_provenance`](https://gdrplatform.github.io/gDRstyle/reference/refit_from_provenance.md)
can call it again as `fn(se, ...)`. Functions passed among the arguments
(e.g. a score function) are stored as function records.

`record_processing_step()` is for a step whose result cannot be
recomputed from the object it returns, e.g. growth rates computed from
an assay that the output no longer carries. It stores the functions and
settings of the step under `provenance$steps[[step]]`.

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[["single-agent"]]
se <- record_processing_step(
  se, "growth_rates",
  settings = list(periods = list(early = c(0, 48))),
  fns = list(rate_fn = function(dt) 0)
)
names(get_SE_provenance(se)$steps)
#> [1] "growth_rates"
```
