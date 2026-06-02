# Define matrix grid positions

Define matrix grid positions

## Usage

``` r
define_matrix_grid_positions(conc1, conc2)
```

## Arguments

- conc1:

  drug_1 concentration

- conc2:

  drug_2 concentration

## Value

list with axis grid positions

## Details

`drug_1` is diluted along the rows as the y-axis and `drug_2` is diluted
along the columns and will be the x-axis.

## Examples

``` r
cl_name <- "cellline_BC"
drug1_name <- "drug_001"
drug2_name <- "drug_026"

se <- get_synthetic_data("combo_matrix_small")[["combination"]]
dt_average <- convert_se_assay_to_dt(se, "Averaged")[normalization_type == "GR"]

ls_axes <- define_matrix_grid_positions(
   dt_average[["Concentration"]], dt_average[["Concentration_2"]])
```
