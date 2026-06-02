# Check whether or not an assay exists in a SummarizedExperiment object.

Check for the presence of an assay in a SummarizedExperiment object.

## Usage

``` r
validate_se_assay_name(se, name)
```

## Arguments

- se:

  A SummarizedExperiment object.

- name:

  String of name of the assay to validate.

## Value

`NULL` invisibly if the assay name is valid. Throws an error if the
assay is not valid.

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
validate_se_assay_name(se, "RawTreated")
```
