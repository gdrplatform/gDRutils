# Apply a function to every element of a bumpy matrix.

Apply a user-specified function to every element of a bumpy matrix.

## Usage

``` r
apply_bumpy_function(
  se,
  FUN,
  req_assay_name,
  out_assay_name,
  parallelize = FALSE,
  ...
)
```

## Arguments

- se:

  A `SummarizedExperiment` object with bumpy matrices.

- FUN:

  A function that will be applied to each element of the matrix in assay
  `req_assay_name`. Output of the function must return a data.table.

- req_assay_name:

  String of the assay name in the `se` that the `FUN` will act on.

- out_assay_name:

  String of the assay name that will contain the results of the applied
  function.

- parallelize:

  Logical indicating whether or not to parallelize the computation.

- ...:

  Additional args to be passed to teh `FUN`.

## Value

The original `se` object with a new assay, `out_assay_name`.

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
FUN <- function(x) {
  data.table::data.table(Concentration = x$Concentration, CorrectedReadout = x$CorrectedReadout)
}
apply_bumpy_function(
  se,
  FUN = FUN,
  req_assay_name = "RawTreated",
  out_assay_name = "CorrectedReadout"
)
#> class: SummarizedExperiment 
#> dim: 10 10 
#> metadata(5): identifiers experiment_metadata Keys fit_parameters
#>   .internal
#> assays(6): RawTreated Controls ... Metrics CorrectedReadout
#> rownames(10): G00002_drug_002_moa_A_72 G00003_drug_003_moa_A_72 ...
#>   G00010_drug_010_moa_A_72 G00011_drug_011_moa_B_72
#> rowData names(4): Gnumber DrugName drug_moa Duration
#> colnames(10): CL00011_cellline_BA_tissue_x_26
#>   CL00012_cellline_CA_tissue_x_30 ... CL00019_cellline_JB_tissue_z_58
#>   CL00020_cellline_KB_tissue_z_62
#> colData names(4): clid CellLineName Tissue ReferenceDivisionTime
```
