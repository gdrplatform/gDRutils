# df_to_bm_assay

Convert data.table with dose-response data into a BumpyMatrix assay.

## Usage

``` r
df_to_bm_assay(data, discard_keys = NULL, precomputed_metadata = NULL)
```

## Arguments

- data:

  data.table with drug-response data

- discard_keys:

  a vector of keys that should be discarded

- precomputed_metadata:

  optional named list from a prior
  [`split_SE_components`](https://gdrplatform.github.io/gDRstyle/reference/split_SE_components.md)
  call to avoid recomputation

## Value

BumpyMatrix object

## Details

The 'assay' is simply a BumpyMatrix object with rownames being the
treatment ids, colnames being the ids of the cell lines and values with
dose-response data for given cell lines under given conditions.

## Examples

``` r
df_to_bm_assay(data.table::data.table(Gnumber = 2, clid = "A"))
#> 1 x 1 BumpyDataFrameMatrix
#> rownames: 1 
#> colnames: 1 
#> preview [1,1]:
#>   DataFrame with 0 rows and 0 columns
```
