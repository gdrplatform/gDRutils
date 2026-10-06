# Convert a SummarizedExperiment assay to a long data.table

Convert an assay within a SummarizedExperiment object to a long
data.table.

## Usage

``` r
convert_se_assay_to_dt(
  se,
  assay_name,
  include_metadata = TRUE,
  retain_nested_rownames = FALSE,
  wide_structure = FALSE,
  unify_metadata = FALSE,
  drop_masked = TRUE,
  merge_additional_variables = FALSE
)
```

## Arguments

- se:

  A SummarizedExperiment object holding raw and/or processed
  dose-response data in its assays.

- assay_name:

  String of name of the assay to transform within the `se`.

- include_metadata:

  Boolean indicating whether or not to include `rowData(se)` and
  `colData(se)` in the returned data.table. Defaults to `TRUE`.

- retain_nested_rownames:

  Boolean indicating whether or not to retain the rownames nested within
  a `BumpyMatrix` assay. Defaults to `FALSE`. If the `assay_name` is not
  of the `BumpyMatrix` class, this argument's value is ignored. If
  `TRUE`, the resulting column in the data.table will be named as
  `"<assay_name>_rownames"`.

- wide_structure:

  Boolean indicating whether or not to transform data.table into wide
  format. `wide_structure = TRUE` requires
  `retain_nested_rownames = TRUE`.

- unify_metadata:

  Boolean indicating whether to unify DrugName and CellLineName in cases
  where DrugNames and CellLineNames are shared by more than one Gnumber
  and/or clid within the experiment.

- drop_masked:

  Boolean indicating whether to drop masked values; TRUE by default.

- merge_additional_variables:

  Boolean indicating whether to merge additional variables identified by
  `get_additional_variables` into the `DrugName` column. Defaults to
  `FALSE`.

## Value

data.table representation of the data in `assay_name`.

## Details

NOTE: to extract information about 'Control' data, simply call the
function with the name of the assay holding data on controls. To extract
the reference data in to same format as 'Averaged' use
`convert_se_ref_assay_to_dt`.

## See also

flatten

## Examples

``` r
mae <- get_synthetic_data("finalMAE_small.qs2")
se <- mae[[1]]
convert_se_assay_to_dt(se, "Metrics")
#>                           rId                             cId
#>                        <char>                          <char>
#>   1: G00002_drug_002_moa_A_72 CL00011_cellline_BA_tissue_x_26
#>   2: G00002_drug_002_moa_A_72 CL00011_cellline_BA_tissue_x_26
#>   3: G00003_drug_003_moa_A_72 CL00011_cellline_BA_tissue_x_26
#>   4: G00003_drug_003_moa_A_72 CL00011_cellline_BA_tissue_x_26
#>   5: G00004_drug_004_moa_A_72 CL00011_cellline_BA_tissue_x_26
#>  ---                                                         
#> 196: G00009_drug_009_moa_A_72 CL00020_cellline_KB_tissue_z_62
#> 197: G00010_drug_010_moa_A_72 CL00020_cellline_KB_tissue_z_62
#> 198: G00010_drug_010_moa_A_72 CL00020_cellline_KB_tissue_z_62
#> 199: G00011_drug_011_moa_B_72 CL00020_cellline_KB_tissue_z_62
#> 200: G00011_drug_011_moa_B_72 CL00020_cellline_KB_tissue_z_62
#>      normalization_type     x_mean     x_AOC x_AOC_range       x_max   x_sd_avg
#>                  <char>      <num>     <num>       <num>       <num>      <num>
#>   1:                 GR  0.4729298 0.5270702   0.6048472  0.33660000 0.03203268
#>   2:                 RV  0.4466890 0.5533110   0.6277120  0.32793333 0.02424423
#>   3:                 GR  0.3932199 0.6067801   0.7027714  0.25136667 0.03075240
#>   4:                 RV  0.3939872 0.6060128   0.6954994  0.27413333 0.02195104
#>   5:                 GR  0.8393271 0.1606729   0.1905009  0.76383333 0.02000573
#>  ---                                                                           
#> 196:                 RV  0.5045079 0.4954921   0.5777711  0.20156667 0.02298094
#> 197:                 GR  0.1991504 0.8008496   0.8915995 -0.78626667 0.06333498
#> 198:                 RV  0.5769999 0.4230001   0.4723104  0.07513333 0.03151223
#> 199:                 GR -0.1148134 1.1148134   1.3059605 -0.76026667 0.07007818
#> 200:                 RV  0.4117805 0.5882195   0.6900581  0.08603333 0.03335063
#>      N_conc maxlog10Concentration        ec50        xc50        h        r2
#>       <int>                 <num>       <num>       <num>    <num>     <num>
#>   1:      9                     1 0.004512560 0.008956223 1.911901 0.9884406
#>   2:      9                     1 0.003795893 0.007253171 1.827031 0.9924695
#>   3:      9                     1 0.004714084 0.006591468 2.270093 0.9905292
#>   4:      9                     1 0.003962986 0.005655581 2.349599 0.9958511
#>   5:      9                     1 0.017152797         Inf 3.143707 0.9914280
#>  ---                                                                        
#> 196:      9                     1 0.035484367 0.045034344 1.998935 0.9961207
#> 197:      9                     1 0.146270629 0.096729600 2.204685 0.9985705
#> 198:      9                     1 0.137152296 0.150209032 2.245643 0.9988156
#> 199:      9                     1 0.026434030 0.016239553 1.854308 0.9987249
#> 200:      9                     1 0.024485958 0.027440547 1.884916 0.9986615
#>               rss      p_value   x_0       x_inf               fit_type
#>             <num>        <num> <num>       <num>                 <char>
#>   1: 0.0045303452 1.660656e-07     1  0.36516781 DRC3pHillFitModelFixS0
#>   2: 0.0027788222 3.705769e-08     1  0.34682665 DRC3pHillFitModelFixS0
#>   3: 0.0056700663 8.267143e-08     1  0.26639679 DRC3pHillFitModelFixS0
#>   4: 0.0021823404 4.600131e-09     1  0.28319787 DRC3pHillFitModelFixS0
#>   5: 0.0007996589 5.831599e-08     1  0.76697772 DRC3pHillFitModelFixS0
#>  ---                                                                   
#> 196: 0.0042739653 3.636093e-09     1  0.18949638 DRC3pHillFitModelFixS0
#> 197: 0.0081610637 1.104466e-10     1 -0.74430401 DRC3pHillFitModelFixS0
#> 198: 0.0018424204 5.718652e-11     1  0.09235428 DRC3pHillFitModelFixS0
#> 199: 0.0064842414 7.403286e-11     1 -0.73401896 DRC3pHillFitModelFixS0
#> 200: 0.0018446214 8.773134e-11     1  0.09662172 DRC3pHillFitModelFixS0
#>      fit_source Gnumber DrugName drug_moa Duration    clid CellLineName
#>          <char>  <char>   <char>   <char>    <num>  <char>       <char>
#>   1:        gDR  G00002 drug_002    moa_A       72 CL00011  cellline_BA
#>   2:        gDR  G00002 drug_002    moa_A       72 CL00011  cellline_BA
#>   3:        gDR  G00003 drug_003    moa_A       72 CL00011  cellline_BA
#>   4:        gDR  G00003 drug_003    moa_A       72 CL00011  cellline_BA
#>   5:        gDR  G00004 drug_004    moa_A       72 CL00011  cellline_BA
#>  ---                                                                   
#> 196:        gDR  G00009 drug_009    moa_A       72 CL00020  cellline_KB
#> 197:        gDR  G00010 drug_010    moa_A       72 CL00020  cellline_KB
#> 198:        gDR  G00010 drug_010    moa_A       72 CL00020  cellline_KB
#> 199:        gDR  G00011 drug_011    moa_B       72 CL00020  cellline_KB
#> 200:        gDR  G00011 drug_011    moa_B       72 CL00020  cellline_KB
#>        Tissue ReferenceDivisionTime
#>        <char>                 <num>
#>   1: tissue_x                    26
#>   2: tissue_x                    26
#>   3: tissue_x                    26
#>   4: tissue_x                    26
#>   5: tissue_x                    26
#>  ---                               
#> 196: tissue_z                    62
#> 197: tissue_z                    62
#> 198: tissue_z                    62
#> 199: tissue_z                    62
#> 200: tissue_z                    62
```
