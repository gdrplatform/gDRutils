# Cap infinity values (Inf, -Inf) in the assay data

Cap infinity values (Inf, -Inf) in the assay data

## Usage

``` r
cap_assay_infinities(
  conc_assay_dt,
  assay_dt,
  experiment_name,
  col = "xc50",
  capping_fold = 5,
  additional_group_cols = NULL
)
```

## Arguments

- conc_assay_dt:

  assay data in data.table format with Concentration data

- assay_dt:

  assay data in data.table format with infinity values to be capped

- experiment_name:

  string with the name of the experiment

- col:

  string with column name to be capped in assay_dt ("xc50" by default)

- capping_fold:

  number for min and max concentration values final formulas are min /
  capping_fold and max \* capping_fold

- additional_group_cols:

  character vector of column names used to identify unique observations

  - for single-agent experiment additional to the combination of
    `DrugName` and `CellLineName`

  - for combination experiment additional to the combination of
    `DrugName`, `DrugName_2` and `CellLineName`

## Value

data.table without -Inf / Inf values

## Examples

``` r
# single-agent data
sdata <- get_synthetic_data("finalMAE_small.qs2")
smetrics_data <- convert_se_assay_to_dt(sdata[[get_supported_experiments("sa")]], "Metrics")
saveraged_data <- convert_se_assay_to_dt(sdata[[get_supported_experiments("sa")]], "Averaged")
smetrics_data_capped <- cap_assay_infinities(saveraged_data,
                                             smetrics_data,
                                             experiment_name = "single-agent")

# combination data
cdata <- get_synthetic_data("finalMAE_combo_matrix_small.qs2")
scaveraged_data <- convert_se_assay_to_dt(cdata[[get_supported_experiments("combo")]], "Averaged")
scmetrics_data <- convert_se_assay_to_dt(cdata[[get_supported_experiments("combo")]], "Metrics")
scmetrics_data_capped <- cap_assay_infinities(scaveraged_data,
                                              scmetrics_data,
                                              experiment_name = "combination")
```
