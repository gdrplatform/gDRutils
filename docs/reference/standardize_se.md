# Standardize SE by switching from custom identifiers into gDR-default

Standardize SE by switching from custom identifiers into gDR-default

## Usage

``` r
standardize_se(se, use_default = TRUE)
```

## Arguments

- se:

  a SummarizedExperiment object with drug-response data generate by gDR
  pipeline

- use_default:

  boolean indicating whether or not to use default identifiers for
  standardization

## Value

se a SummarizedExperiment with default gDR identifiers

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
S4Vectors::metadata(se)$identifiers$drug <- "druug"
standardize_se(se)
#> Warning: overwriting existing metadata entry: 'identifiers'
#> class: SummarizedExperiment 
#> dim: 10 10 
#> metadata(5): identifiers experiment_metadata Keys fit_parameters
#>   .internal
#> assays(5): RawTreated Controls Normalized Averaged Metrics
#> rownames(10): G00002_drug_002_moa_A_72 G00003_drug_003_moa_A_72 ...
#>   G00010_drug_010_moa_A_72 G00011_drug_011_moa_B_72
#> rowData names(4): Gnumber DrugName drug_moa Duration
#> colnames(10): CL00011_cellline_BA_tissue_x_26
#>   CL00012_cellline_CA_tissue_x_30 ... CL00019_cellline_JB_tissue_z_58
#>   CL00020_cellline_KB_tissue_z_62
#> colData names(4): clid CellLineName Tissue ReferenceDivisionTime
```
