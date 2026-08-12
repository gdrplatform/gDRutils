# Merge assay data

Merge assay data

## Usage

``` r
merge_assay(
  SElist,
  assay_name,
  additional_col_name = "data_source",
  discard_keys = NULL
)
```

## Arguments

- SElist:

  named list of Summarized Experiments

- assay_name:

  name of the assay that should be extracted and merged

- additional_col_name:

  string of column name that will be added to assay data for the
  distinction of possible duplicated metrics that can arise from
  multiple projects

- discard_keys:

  character vector of string that will be discarded during creating
  BumpyMatrix object

## Value

BumpyMatrix or list with data.table + BumpyMatrix

## Examples

``` r
mae <- get_synthetic_data("finalMAE_combo_2dose_nonoise.qs2")

listSE <- list(
  combo1 = mae[[1]],
  sa = mae[[2]]
)
merge_assay(listSE, "Normalized")
#> $DT
#>       data_source Concentration Concentration_2 normalization_type     x
#>            <char>         <num>           <num>             <char> <num>
#>    1:      combo1   0.001000000             0.0                 RV 0.957
#>    2:      combo1   0.001000000             0.2                 RV 0.957
#>    3:      combo1   0.001000000             1.0                 RV 0.957
#>    4:      combo1   0.003162278             0.0                 RV 0.712
#>    5:      combo1   0.003162278             0.2                 RV 0.712
#>   ---                                                                   
#> 2804:          sa   1.000000000              NA                 GR 1.000
#> 2805:          sa   1.000000000              NA                 GR 1.000
#> 2806:          sa   1.000000000              NA                 GR 1.000
#> 2807:          sa   1.000000000              NA                 GR 1.000
#> 2808:          sa   1.000000000              NA                 GR 1.000
#>       Gnumber DrugName drug_moa Gnumber_2 DrugName_2 drug_moa_2 Duration
#>        <char>   <char>   <char>    <char>     <char>     <char>    <num>
#>    1:  G00002 drug_002    moa_A    G00026   drug_026      moa_E       72
#>    2:  G00002 drug_002    moa_A    G00026   drug_026      moa_E       72
#>    3:  G00002 drug_002    moa_A    G00026   drug_026      moa_E       72
#>    4:  G00002 drug_002    moa_A    G00026   drug_026      moa_E       72
#>    5:  G00002 drug_002    moa_A    G00026   drug_026      moa_E       72
#>   ---                                                                   
#> 2804:  G00026 drug_026    moa_E      <NA>       <NA>       <NA>       72
#> 2805:  G00026 drug_026    moa_E      <NA>       <NA>       <NA>       72
#> 2806:  G00026 drug_026    moa_E      <NA>       <NA>       <NA>       72
#> 2807:  G00026 drug_026    moa_E      <NA>       <NA>       <NA>       72
#> 2808:  G00026 drug_026    moa_E      <NA>       <NA>       <NA>       72
#>          clid CellLineName   Tissue ReferenceDivisionTime
#>        <char>       <char>   <char>                 <num>
#>    1: CL00011  cellline_BA tissue_x                    26
#>    2: CL00011  cellline_BA tissue_x                    26
#>    3: CL00011  cellline_BA tissue_x                    26
#>    4: CL00011  cellline_BA tissue_x                    26
#>    5: CL00011  cellline_BA tissue_x                    26
#>   ---                                                    
#> 2804: CL00013  cellline_DA tissue_x                    34
#> 2805: CL00013  cellline_DA tissue_x                    34
#> 2806: CL00013  cellline_DA tissue_x                    34
#> 2807: CL00013  cellline_DA tissue_x                    34
#> 2808: CL00013  cellline_DA tissue_x                    34
#> 
#> $BM
#> 13 x 3 BumpyDataFrameMatrix
#> rownames: 1 2 ... 12 13 
#> colnames: 1 2 3 
#> preview [1,1]:
#>   DataFrame with 54 rows and 4 columns
#>               x Concentration data_source normalization_type
#>       <numeric>     <numeric> <character>        <character>
#>   1      0.9570    0.00100000      combo1                 RV
#>   2      0.7120    0.00316228      combo1                 RV
#>   3      0.4119    0.01000000      combo1                 RV
#>   4      0.3501    0.03162278      combo1                 RV
#>   5      0.3317    0.10000000      combo1                 RV
#>   ...       ...           ...         ...                ...
#>   50     0.3427      0.100000      combo1                 GR
#>   51     0.2086      0.316228      combo1                 GR
#>   52     0.1019      1.000000      combo1                 GR
#>   53     0.0908      3.162278      combo1                 GR
#>   54     0.0902     10.000000      combo1                 GR
#> 
```
