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
#>                           rId                             cId     x_mean
#>                        <char>                          <char>      <num>
#>   1: G00002_drug_002_moa_A_72 CL00011_cellline_BA_tissue_x_26  0.4466890
#>   2: G00002_drug_002_moa_A_72 CL00011_cellline_BA_tissue_x_26  0.4729298
#>   3: G00003_drug_003_moa_A_72 CL00011_cellline_BA_tissue_x_26  0.3939872
#>   4: G00003_drug_003_moa_A_72 CL00011_cellline_BA_tissue_x_26  0.3932199
#>   5: G00004_drug_004_moa_A_72 CL00011_cellline_BA_tissue_x_26  0.7989283
#>  ---                                                                    
#> 196: G00009_drug_009_moa_A_72 CL00020_cellline_KB_tissue_z_62  0.0774842
#> 197: G00010_drug_010_moa_A_72 CL00020_cellline_KB_tissue_z_62  0.5769999
#> 198: G00010_drug_010_moa_A_72 CL00020_cellline_KB_tissue_z_62  0.1991504
#> 199: G00011_drug_011_moa_B_72 CL00020_cellline_KB_tissue_z_62  0.4117805
#> 200: G00011_drug_011_moa_B_72 CL00020_cellline_KB_tissue_z_62 -0.1148134
#>          x_AOC x_AOC_range        xc50       x_max        ec50       x_inf
#>          <num>       <num>       <num>       <num>       <num>       <num>
#>   1: 0.5533110   0.6277120 0.007253171  0.32793333 0.003795893  0.34682665
#>   2: 0.5270702   0.6048472 0.008956223  0.33660000 0.004512560  0.36516781
#>   3: 0.6060128   0.6954994 0.005655581  0.27413333 0.003962986  0.28319787
#>   4: 0.6067801   0.7027714 0.006591468  0.25136667 0.004714084  0.26639679
#>   5: 0.2010717   0.2385432         Inf  0.70636667 0.016546380  0.71001484
#>  ---                                                                      
#> 196: 0.9225158   1.0742560 0.026077313 -0.49686667 0.037516853 -0.52393526
#> 197: 0.4230001   0.4723104 0.150209032  0.07513333 0.137152296  0.09235428
#> 198: 0.8008496   0.8915995 0.096729600 -0.78626667 0.146270629 -0.74430401
#> 199: 0.5882195   0.6900581 0.027440547  0.08603333 0.024485958  0.09662172
#> 200: 1.1148134   1.3059605 0.016239553 -0.76026667 0.026434030 -0.73401896
#>        x_0        h        r2      p_value         rss maxlog10Concentration
#>      <num>    <num>     <num>        <num>       <num>                 <num>
#>   1:     1 1.827031 0.9924695 3.705769e-08 0.002778822                     1
#>   2:     1 1.911901 0.9884406 1.660656e-07 0.004530345                     1
#>   3:     1 2.349599 0.9958511 4.600131e-09 0.002182340                     1
#>   4:     1 2.270093 0.9905292 8.267143e-08 0.005670066                     1
#>   5:     1 3.159939 0.9918874 4.809026e-08 0.001162341                     1
#>  ---                                                                        
#> 196:     1 1.970723 0.9955609 5.828304e-09 0.017400921                     1
#> 197:     1 2.245643 0.9988156 5.718652e-11 0.001842420                     1
#> 198:     1 2.204685 0.9985705 1.104466e-10 0.008161064                     1
#> 199:     1 1.884916 0.9986615 8.773134e-11 0.001844621                     1
#> 200:     1 1.854308 0.9987249 7.403286e-11 0.006484241                     1
#>      N_conc   x_sd_avg               fit_type normalization_type fit_source
#>       <int>      <num>                 <char>             <char>     <char>
#>   1:      9 0.02424423 DRC3pHillFitModelFixS0                 RV        gDR
#>   2:      9 0.03203268 DRC3pHillFitModelFixS0                 GR        gDR
#>   3:      9 0.02195104 DRC3pHillFitModelFixS0                 RV        gDR
#>   4:      9 0.03075240 DRC3pHillFitModelFixS0                 GR        gDR
#>   5:      9 0.02315049 DRC3pHillFitModelFixS0                 RV        gDR
#>  ---                                                                       
#> 196:      9 0.04602459 DRC3pHillFitModelFixS0                 GR        gDR
#> 197:      9 0.03151223 DRC3pHillFitModelFixS0                 RV        gDR
#> 198:      9 0.06333498 DRC3pHillFitModelFixS0                 GR        gDR
#> 199:      9 0.03335063 DRC3pHillFitModelFixS0                 RV        gDR
#> 200:      9 0.07007818 DRC3pHillFitModelFixS0                 GR        gDR
#>      Gnumber DrugName drug_moa Duration    clid CellLineName   Tissue
#>       <char>   <char>   <char>    <num>  <char>       <char>   <char>
#>   1:  G00002 drug_002    moa_A       72 CL00011  cellline_BA tissue_x
#>   2:  G00002 drug_002    moa_A       72 CL00011  cellline_BA tissue_x
#>   3:  G00003 drug_003    moa_A       72 CL00011  cellline_BA tissue_x
#>   4:  G00003 drug_003    moa_A       72 CL00011  cellline_BA tissue_x
#>   5:  G00004 drug_004    moa_A       72 CL00011  cellline_BA tissue_x
#>  ---                                                                 
#> 196:  G00009 drug_009    moa_A       72 CL00020  cellline_KB tissue_z
#> 197:  G00010 drug_010    moa_A       72 CL00020  cellline_KB tissue_z
#> 198:  G00010 drug_010    moa_A       72 CL00020  cellline_KB tissue_z
#> 199:  G00011 drug_011    moa_B       72 CL00020  cellline_KB tissue_z
#> 200:  G00011 drug_011    moa_B       72 CL00020  cellline_KB tissue_z
#>      ReferenceDivisionTime
#>                      <num>
#>   1:                    26
#>   2:                    26
#>   3:                    26
#>   4:                    26
#>   5:                    26
#>  ---                      
#> 196:                    62
#> 197:                    62
#> 198:                    62
#> 199:                    62
#> 200:                    62
```
