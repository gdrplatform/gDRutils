# mcolData

get colData of all experiments

## Usage

``` r
mcolData(mae)
```

## Arguments

- mae:

  MultiAssayExperiment object

## Value

data.table with all-experiments colData

## Author

Arkadiusz Gladki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
mcolData(mae)
#>        clid CellLineName   Tissue ReferenceDivisionTime
#>      <char>       <char>   <char>                 <num>
#>  1: CL00011  cellline_BA tissue_x                    26
#>  2: CL00012  cellline_CA tissue_x                    30
#>  3: CL00013  cellline_DA tissue_x                    34
#>  4: CL00014  cellline_EA tissue_x                    38
#>  5: CL00015  cellline_FA tissue_x                    42
#>  6: CL00016  cellline_GB tissue_y                    46
#>  7: CL00017  cellline_HB tissue_y                    50
#>  8: CL00018  cellline_IB tissue_y                    54
#>  9: CL00019  cellline_JB tissue_z                    58
#> 10: CL00020  cellline_KB tissue_z                    62
```
