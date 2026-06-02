# Create JSON document.

Convert a MultiAssayExperiment object to a JSON document.

## Usage

``` r
convert_mae_to_json(mae, with_experiments = TRUE)
```

## Arguments

- mae:

  SummarizedExperiment object.

- with_experiments:

  logical convert experiment metadata as well?

## Value

String representation of a JSON document.

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
convert_mae_to_json(mae)
#> $mae
#> {"experiment_names":["single-agent"]} 
#> 
#> $se
#> $se$`single-agent`
#> {"drug":["G00002","G00003","G00004","G00005","G00006","G00007","G00008","G00009","G00010","G00011"],"drug_name":["drug_002","drug_003","drug_004","drug_005","drug_006","drug_007","drug_008","drug_009","drug_010","drug_011"],"drug_moa":["moa_A","moa_A","moa_A","moa_A","moa_A","moa_A","moa_A","moa_A","moa_A","moa_B"],"duration":[72,72,72,72,72,72,72,72,72,72],"misc_rowdata":{},"cellline":["CL00011","CL00012","CL00013","CL00014","CL00015","CL00016","CL00017","CL00018","CL00019","CL00020"],"cellline_name":["cellline_BA","cellline_CA","cellline_DA","cellline_EA","cellline_FA","cellline_GB","cellline_HB","cellline_IB","cellline_JB","cellline_KB"],"cellline_tissue":["tissue_x","tissue_x","tissue_x","tissue_x","tissue_x","tissue_y","tissue_y","tissue_y","tissue_z","tissue_z"],"cellline_ref_div_time":[26,30,34,38,42,46,50,54,58,62],"misc_coldata":{}} 
#> 
#> 
convert_mae_to_json(mae, with_experiments = FALSE)
#> $mae
#> {"experiment_names":["single-agent"]} 
#> 
```
