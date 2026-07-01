# Conditional lapply with optional batch processing.

Conditional lapply with optional batch processing.

## Usage

``` r
loop(
  x,
  FUN,
  parallelize = as.logical(Sys.getenv("GDR_PARALLELIZE", "FALSE")),
  use_batch = as.logical(Sys.getenv("GDR_USE_BATCH", "FALSE")),
  temp_dir = Sys.getenv("GDR_TEMP_DIR", tempdir()),
  batch_size = as.numeric(Sys.getenv("GDR_BATCH_SIZE", 100)),
  ...
)
```

## Arguments

- x:

  Vector (atomic or list) or an expression object. Other objects
  (including classed objects) will be coerced by
  [as.list](https://rdrr.io/r/base/list.html)

- FUN:

  A user-defined function to apply to each element of `x`.

- parallelize:

  Logical indicating whether or not to parallelize the computation.
  Defaults to `as.logical(Sys.getenv("GDR_PARALLELIZE", "FALSE"))`.

- use_batch:

  Logical indicating whether to use batch processing to save
  intermediate results. Defaults to `FALSE`.

- temp_dir:

  Character string specifying the directory where batch results are
  saved. Defaults to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- batch_size:

  Integer specifying the number of elements to process in each batch
  during batch mode. Defaults to `100`.

- ...:

  Optional arguments passed to
  [lapply](https://rdrr.io/r/base/lapply.html).

## Value

List containing output of `FUN` applied to every element in `x`. When
batch processing is enabled, results are saved incrementally and merged
at the end of processing.

## Details

The function operates in two modes:

1.  Regular mode: Directly applies `FUN` to the elements using `lapply`.

2.  Batch mode: Saves results in batches to disk, allowing computation
    to resume from the last saved step. Batch mode is activated by
    setting `use_batch` to `TRUE`.

## Examples

``` r
# Regular processing
loop(list(1, 2, 3), function(x) x^2, parallelize = FALSE, use_batch = FALSE)

# Batch processing
loop(1:10, function(x) x^2, parallelize = TRUE, use_batch = TRUE)
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 4
#> 
#> [[3]]
#> [1] 9
#> 
#> [[4]]
#> [1] 16
#> 
#> [[5]]
#> [1] 25
#> 
#> [[6]]
#> [1] 36
#> 
#> [[7]]
#> [1] 49
#> 
#> [[8]]
#> [1] 64
#> 
#> [[9]]
#> [1] 81
#> 
#> [[10]]
#> [1] 100
#> 
```
