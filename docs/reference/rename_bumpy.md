# Rename BumpyMatrix

Rename BumpyMatrix

## Usage

``` r
rename_bumpy(bumpy, mapping_vector)
```

## Arguments

- bumpy:

  a BumpyMatrix object

- mapping_vector:

  a named vector for mapping old and new values. The names of the
  character vector indicate the source names, and the corresponding
  values the destination names.

## Value

a renamed BumpyMatrix object

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
assay <- SummarizedExperiment::assay(se)
rename_bumpy(assay, c("Concentration" = "conc"))
#> 10 x 10 BumpyDataFrameMatrix
#> rownames: G00002_drug_002_moa_A_72 G00003_drug_003_moa_A_72 ... G00010_drug_010_moa_A_72 G00011_drug_011_moa_B_72 
#> colnames: CL00011_cellline_BA_tissue_x_26 CL00012_cellline_CA_tissue_x_30 ... CL00019_cellline_JB_tissue_z_58 CL00020_cellline_KB_tissue_z_62 
#> preview [1,1]:
#>   DataFrame with 27 rows and 6 columns
#>           Barcode       conc ReadoutValue BackgroundValue record_id
#>       <character>  <numeric>    <numeric>       <numeric> <integer>
#>   1       plate_1 0.00100000         93.5               0       601
#>   2       plate_2 0.00100000         92.0               0       701
#>   3       plate_3 0.00100000         91.2               0       801
#>   4       plate_1 0.00316228         74.8               0       901
#>   5       plate_2 0.00316228         74.3               0      1001
#>   ...         ...        ...          ...             ...       ...
#>   23      plate_2    3.16228         38.6               0      2801
#>   24      plate_3    3.16228         33.1               0      2901
#>   25      plate_1   10.00000         31.9               0      3001
#>   26      plate_2   10.00000         30.8               0      3101
#>   27      plate_3   10.00000         35.4               0      3201
#>       CorrectedReadout
#>              <numeric>
#>   1               93.5
#>   2               92.0
#>   3               91.2
#>   4               74.8
#>   5               74.3
#>   ...              ...
#>   23              38.6
#>   24              33.1
#>   25              31.9
#>   26              30.8
#>   27              35.4
```
