# Merge multiple MultiAssayExperiment objects

Merge multiple MultiAssayExperiment objects

## Usage

``` r
merge_MAE(
  MAElist,
  additional_col_name = "data_source",
  discard_keys = c("normalization_type", "fit_source", "record_id", "isDay0", "swap_sa",
    "control_type", "iso_level", "conc_1", "conc_2"),
  title = NULL,
  description = NULL,
  source_name = NULL,
  source_id = NULL
)
```

## Arguments

- MAElist:

  Named list of MultiAssayExperiment objects.

- additional_col_name:

  String with the name of the column that will be added to assay data
  for the distinction of possible duplicated metrics that can arise from
  multiple projects.

- discard_keys:

  Character vector of strings that will be discarded during creating
  BumpyMatrix object.

- title:

  String specifying the final DataSetDB title. If NULL, auto-generates.

- description:

  String specifying the final DataSetDB description. If NULL,
  auto-generates.

- source_name:

  String specifying the standard DSDB source name. If NULL, auto-detects
  or uses "merged_analysis".

- source_id:

  String specifying the unique DSDB source ID. If NULL, uses
  "merged_dataset".

## Value

Merged MultiAssayExperiment object.

## Examples

``` r
mae1 <- get_synthetic_data("finalMAE_combo_2dose_nonoise.qs2")
mae2 <- get_synthetic_data("finalMAE_combo_2dose_nonoise.qs2")
merge_MAE(list(mae1 = mae1, mae2 = mae2), title = "Test", description = "Test MAE")
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] combination: SummarizedExperiment with 3 rows and 3 columns
#>  [2] single-agent: SummarizedExperiment with 4 rows and 3 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
