# mrowData

get rowData of all experiments

## Usage

``` r
mrowData(mae)
```

## Arguments

- mae:

  MultiAssayExperiment object

## Value

data.table with all-experiments rowData

## Author

Arkadiusz Gladki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
mrowData(mae)
#>     Gnumber DrugName drug_moa Duration
#>      <char>   <char>   <char>    <num>
#>  1:  G00002 drug_002    moa_A       72
#>  2:  G00003 drug_003    moa_A       72
#>  3:  G00004 drug_004    moa_A       72
#>  4:  G00005 drug_005    moa_A       72
#>  5:  G00006 drug_006    moa_A       72
#>  6:  G00007 drug_007    moa_A       72
#>  7:  G00008 drug_008    moa_A       72
#>  8:  G00009 drug_009    moa_A       72
#>  9:  G00010 drug_010    moa_A       72
#> 10:  G00011 drug_011    moa_B       72
```
