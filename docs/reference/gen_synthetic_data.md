# gen_synthetic_data

Function for generating local synthetic data used for unit tests in
modules

## Usage

``` r
gen_synthetic_data(m = 1, n = 5)
```

## Arguments

- m:

  number of drugs

- n:

  number of records

## Value

list with drugs, cell_lines, raw_data and assay_data

## Examples

``` r
gen_synthetic_data()
#> $drug_names
#> [1] "drug_001"
#> 
#> $cell_names
#> [1] "cellline_BA" "cellline_CA" "cellline_DA" "cellline_EA" "cellline_FA"
#> 
#> $dt
#>    Drug Name Drug MOA Cell Line Name   Tissue GR_AOC GR Inf  GR 0 GEC50  h GR
#>       <char>   <char>         <char>   <char>  <num>  <num> <num> <num> <num>
#> 1:  drug_001    moa_B    cellline_BA tissue_x    0.1    0.1   0.1   0.1   0.1
#> 2:  drug_001    moa_B    cellline_CA tissue_x    0.7    0.7   0.7   0.7   0.7
#> 3:  drug_001    moa_B    cellline_DA tissue_x    1.3    1.3   1.3   1.3   1.3
#> 4:  drug_001    moa_B    cellline_EA tissue_x    1.9    1.9   1.9   1.9   1.9
#> 5:  drug_001    moa_B    cellline_FA tissue_y    2.5    2.5   2.5   2.5   2.5
#>    E Inf    E0  EC50  h RV  GR50  IC50 GR Max E Max GR value Concentration
#>    <num> <num> <num> <num> <num> <num>  <num> <num>    <num>         <num>
#> 1:   0.1   0.1   0.1   0.1   0.1   0.1    0.1   0.1      0.1           0.1
#> 2:   0.7   0.7   0.7   0.7   0.7   0.7    0.7   0.7      0.7           0.7
#> 3:   1.3   1.3   1.3   1.3   1.3   1.3    1.3   1.3      1.3           1.3
#> 4:   1.9   1.9   1.9   1.9   1.9   1.9    1.9   1.9      1.9           1.9
#> 5:   2.5   2.5   2.5   2.5   2.5   2.5    2.5   2.5      2.5           2.5
#> 
```
