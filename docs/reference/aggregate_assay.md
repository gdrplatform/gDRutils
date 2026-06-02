# Aggregate a `BumpyMatrix` assay by a given aggregation function.

Aggregation can only be performed on nested variables.

## Usage

``` r
aggregate_assay(asy, by, FUN)
```

## Arguments

- asy:

  A `BumpyMatrix` object.

- by:

  Character vector of the nested fields to aggregate by.

- FUN:

  A function to use to aggregate the data.

## Value

A `BumpyMatrix` object aggregated by `FUN`.

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
assay <- SummarizedExperiment::assay(se)
#> Loading required namespace: BumpyMatrix
aggregate_assay(assay, FUN = mean, by = c("Barcode"))
#> 10 x 10 BumpyDataFrameMatrix
#> rownames: G00002_drug_002_moa_A_72 G00003_drug_003_moa_A_72 ... G00010_drug_010_moa_A_72 G00011_drug_011_moa_B_72 
#> colnames: CL00011_cellline_BA_tissue_x_26 CL00012_cellline_CA_tissue_x_30 ... CL00019_cellline_JB_tissue_z_58 CL00020_cellline_KB_tissue_z_62 
#> preview [1,1]:
#>   DataFrame with 3 rows and 6 columns
#>         Barcode Concentration ReadoutValue BackgroundValue record_id
#>     <character>     <numeric>    <numeric>       <numeric> <numeric>
#>   1     plate_1       1.62492      46.0333               0      1801
#>   2     plate_2       1.62492      46.4111               0      1901
#>   3     plate_3       1.62492      47.0778               0      2001
#>     CorrectedReadout
#>            <numeric>
#>   1          46.0333
#>   2          46.4111
#>   3          47.0778
```
