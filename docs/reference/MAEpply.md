# Lapply through all the experiments in MultiAssayExperiment object

Lapply through all the experiments in MultiAssayExperiment object

## Usage

``` r
MAEpply(mae, FUN, unify = FALSE, ...)
```

## Arguments

- mae:

  MultiAssayExperiment object

- FUN:

  function that should be applied on each experiment of
  MultiAssayExperiment object

- unify:

  logical indicating if the output should be a unlisted object of unique
  values across all the experiments

- ...:

  Additional args to be passed to teh `FUN`.

## Value

list or vector depends on unify param

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
#> Loading required namespace: MultiAssayExperiment
MAEpply(mae, SummarizedExperiment::assayNames)
#> $`single-agent`
#> [1] "RawTreated" "Controls"   "Normalized" "Averaged"   "Metrics"   
#> 
```
