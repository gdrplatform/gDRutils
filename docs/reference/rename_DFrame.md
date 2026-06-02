# Rename DFrame

Rename DFrame

## Usage

``` r
rename_DFrame(df, mapping_vector)
```

## Arguments

- df:

  a DFrame object

- mapping_vector:

  a named vector for mapping old and new values. The names of the
  character vector indicate the source names, and the corresponding
  values the destination names.

## Value

a renamed DFrame object

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
rename_DFrame(SummarizedExperiment::rowData(mae[[1]]), c("Gnumber" = "Gnumber1"))
#> DataFrame with 10 rows and 4 columns
#>                             Gnumber1    DrugName    drug_moa  Duration
#>                          <character> <character> <character> <numeric>
#> G00002_drug_002_moa_A_72      G00002    drug_002       moa_A        72
#> G00003_drug_003_moa_A_72      G00003    drug_003       moa_A        72
#> G00004_drug_004_moa_A_72      G00004    drug_004       moa_A        72
#> G00005_drug_005_moa_A_72      G00005    drug_005       moa_A        72
#> G00006_drug_006_moa_A_72      G00006    drug_006       moa_A        72
#> G00007_drug_007_moa_A_72      G00007    drug_007       moa_A        72
#> G00008_drug_008_moa_A_72      G00008    drug_008       moa_A        72
#> G00009_drug_009_moa_A_72      G00009    drug_009       moa_A        72
#> G00010_drug_010_moa_A_72      G00010    drug_010       moa_A        72
#> G00011_drug_011_moa_B_72      G00011    drug_011       moa_B        72
```
