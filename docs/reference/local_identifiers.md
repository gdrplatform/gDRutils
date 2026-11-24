# Use a set of identifiers for the duration of a call

`local_identifiers()` makes `ids` the active identifiers until the
calling function returns, on success or on error, then restores the
identifiers that were active before. `with_identifiers()` does the same
around one expression, and `with_SE_identifiers()` uses the identifiers
stored on a `SummarizedExperiment`. Keys absent from `ids` keep their
current value.

## Usage

``` r
local_identifiers(ids, .local_envir = parent.frame())

with_identifiers(ids, code)

with_SE_identifiers(se, code)
```

## Arguments

- ids:

  named list of identifiers, a subset of
  `names(get_default_identifiers())`.

- .local_envir:

  environment whose exit restores the previous identifiers.

- code:

  expression evaluated with the identifiers active.

- se:

  a `SummarizedExperiment`.

## Value

`local_identifiers()`: the identifiers active before, invisibly.
`with_identifiers()` and `with_SE_identifiers()`: the value of `code`.

## Details

Prefer these to
[`set_env_identifier()`](https://gdrplatform.github.io/gDRstyle/reference/identifiers.md),
which changes the identifiers for the rest of the session. Called at the
top level (the console), there is no function to return from, so
`local_identifiers()` warns and the identifiers stay set, as with
[`set_env_identifier()`](https://gdrplatform.github.io/gDRstyle/reference/identifiers.md).

An SE built by gDRcore stores the identifiers it was built with, and the
functions that read an SE
([`convert_se_assay_to_dt()`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_assay_to_dt.md),
[`validate_SE()`](https://gdrplatform.github.io/gDRstyle/reference/validate_SE.md),
[`apply_fit()`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.md),
...) use them in preference to the session's. A key the SE does not
store is taken from the session, as
[`get_SE_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md)
always did. To correct an object's identifiers, use
[`set_SE_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/SE_metadata.md).

## Examples

``` r
f <- function() {
  local_identifiers(list(drug = "DrugID"))
  get_env_identifiers("drug")
}
f()                         # "DrugID"
#> [1] "DrugID"
get_env_identifiers("drug") # unchanged
#> [1] "Gnumber"
with_identifiers(list(cellline = "Clid"), get_env_identifiers("cellline"))
#> [1] "Clid"
```
