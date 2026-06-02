# Prepare dict with min and max concentration for codilution

Prepare dict with min and max concentration for codilution

## Usage

``` r
.prep_cd_conc_cap_dict(
  conc_assay_dt,
  group_cols = as.character(get_env_identifiers(c("drug_name", "drug_name2",
    "cellline_name"), simplify = FALSE))
)
```

## Arguments

- conc_assay_dt:

  assay data in data.table format with Concentration data

- group_cols:

  charvec with grouping column names

## Value

`data.table` with max and min concentration for codilution
