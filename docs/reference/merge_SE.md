# Merge multiple Summarized Experiments

Merge multiple Summarized Experiments

## Usage

``` r
merge_SE(
  SElist,
  additional_col_name = "data_source",
  discard_keys = c("normalization_type", "fit_source", "record_id", "isDay0", "swap_sa",
    "control_type", "iso_level", "conc_1", "conc_2")
)
```

## Arguments

- SElist:

  named list of Summarized Experiments

- additional_col_name:

  string with the name of the column that will be added to assay data
  for the distinction of possible duplicated metrics that can arise from
  multiple projects

- discard_keys:

  character vector of string that will be discarded during creating
  BumpyMatrix object

## Value

merged SummarizedExperiment object

## Examples

``` r
se1 <- get_synthetic_data("finalMAE_small.qs2")[[1]]
merge_SE(list(se1 = se1, se2 = se1))
#> class: SummarizedExperiment 
#> dim: 10 10 
#> metadata(5): experiment_metadata Keys fit_parameters .internal
#>   identifiers
#> assays(5): RawTreated Controls Normalized Averaged Metrics
#> rownames(10): G00002_drug_002_moa_A_72 G00003_drug_003_moa_A_72 ...
#>   G00010_drug_010_moa_A_72 G00011_drug_011_moa_B_72
#> rowData names(4): Gnumber DrugName drug_moa Duration
#> colnames(10): CL00011_cellline_BA_tissue_x_26
#>   CL00012_cellline_CA_tissue_x_30 ... CL00019_cellline_JB_tissue_z_58
#>   CL00020_cellline_KB_tissue_z_62
#> colData names(4): clid CellLineName Tissue ReferenceDivisionTime
```
