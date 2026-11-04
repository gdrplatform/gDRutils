# convert combo assays from SummarizedExperiments to the list of data.tables

convert combo assays from SummarizedExperiments to the list of
data.tables

## Usage

``` r
convert_combo_data_to_dt(
  se,
  c_assays = get_combo_assay_names(),
  normalization_type = c("RV", "GR"),
  prettify = TRUE
)
```

## Arguments

- se:

  `SummarizedExperiment` object with dose-response data

- c_assays:

  charvec of combo assays to be used

- normalization_type:

  charvec of normalization_types expected in the data

- prettify:

  boolean flag indicating whether or not to prettify the colnames of the
  returned data

## Value

list of data.table(s) with combo data

## Author

Arkadiusz Gładki <arkadiusz.gladki@contractors.roche.com>

## Examples

``` r
mae <- get_synthetic_data("finalMAE_combo_matrix_small.qs2")
convert_combo_data_to_dt(mae[[1]])
#> $excess
#>                                                 r Id
#>                                               <char>
#>    1: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>    2: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>    3: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>    4: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>    5: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>   ---                                               
#> 1940: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 1941: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 1942: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 1943: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 1944: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#>                                  c Id Concentration Concentration 2    Smooth
#>                                <char>         <num>           <num>     <num>
#>    1: CL00016_cellline_GB_tissue_y_46          0.00         0.00000 1.0000000
#>    2: CL00016_cellline_GB_tissue_y_46          0.00         0.00100 1.0000000
#>    3: CL00016_cellline_GB_tissue_y_46          0.00         0.00316 1.0000000
#>    4: CL00016_cellline_GB_tissue_y_46          0.00         0.01000 1.0000000
#>    5: CL00016_cellline_GB_tissue_y_46          0.00         0.03160 1.0000000
#>   ---                                                                        
#> 1940: CL00017_cellline_HB_tissue_y_50          3.16         0.03160 0.1737447
#> 1941: CL00017_cellline_HB_tissue_y_50          3.16         0.10000 0.1748643
#> 1942: CL00017_cellline_HB_tissue_y_50          3.16         0.31600 0.1828034
#> 1943: CL00017_cellline_HB_tissue_y_50          3.16         1.00000 0.1918578
#> 1944: CL00017_cellline_HB_tissue_y_50          3.16         3.16000 0.1921170
#>       Normalization Type HSA Excess Bliss Excess Gnumber Drug Name Drug MOA
#>                   <char>      <num>        <num>  <char>    <char>   <char>
#>    1:                 GR         NA           NA  G00004  drug_004    moa_A
#>    2:                 GR         NA           NA  G00004  drug_004    moa_A
#>    3:                 GR         NA           NA  G00004  drug_004    moa_A
#>    4:                 GR         NA           NA  G00004  drug_004    moa_A
#>    5:                 GR         NA           NA  G00004  drug_004    moa_A
#>   ---                                                                      
#> 1940:                 RV 0.06097726   0.06097726  G00006  drug_006    moa_A
#> 1941:                 RV 0.05985762   0.05985762  G00006  drug_006    moa_A
#> 1942:                 RV 0.05191848   0.05191848  G00006  drug_006    moa_A
#> 1943:                 RV 0.04286414   0.04286414  G00006  drug_006    moa_A
#> 1944:                 RV 0.04260493   0.04260493  G00006  drug_006    moa_A
#>       Gnumber 2 Drug Name 2 Drug MOA 2 Duration    Clid Cell Line Name   Tissue
#>          <char>      <char>     <char>    <num>  <char>         <char>   <char>
#>    1:    G00021    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#>    2:    G00021    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#>    3:    G00021    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#>    4:    G00021    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#>    5:    G00021    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#>   ---                                                                          
#> 1940:    G00026    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#> 1941:    G00026    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#> 1942:    G00026    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#> 1943:    G00026    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#> 1944:    G00026    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#>       Reference Division Time
#>                         <num>
#>    1:                      46
#>    2:                      46
#>    3:                      46
#>    4:                      46
#>    5:                      46
#>   ---                        
#> 1940:                      50
#> 1941:                      50
#> 1942:                      50
#> 1943:                      50
#> 1944:                      50
#> 
#> $scores
#>                                               r Id
#>                                             <char>
#>  1: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>  2: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>  3: G00004_drug_004_moa_A_G00026_drug_026_moa_E_72
#>  4: G00004_drug_004_moa_A_G00026_drug_026_moa_E_72
#>  5: G00005_drug_005_moa_A_G00021_drug_021_moa_D_72
#>  6: G00005_drug_005_moa_A_G00021_drug_021_moa_D_72
#>  7: G00005_drug_005_moa_A_G00026_drug_026_moa_E_72
#>  8: G00005_drug_005_moa_A_G00026_drug_026_moa_E_72
#>  9: G00006_drug_006_moa_A_G00021_drug_021_moa_D_72
#> 10: G00006_drug_006_moa_A_G00021_drug_021_moa_D_72
#> 11: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 12: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 13: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#> 14: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#> 15: G00004_drug_004_moa_A_G00026_drug_026_moa_E_72
#> 16: G00004_drug_004_moa_A_G00026_drug_026_moa_E_72
#> 17: G00005_drug_005_moa_A_G00021_drug_021_moa_D_72
#> 18: G00005_drug_005_moa_A_G00021_drug_021_moa_D_72
#> 19: G00005_drug_005_moa_A_G00026_drug_026_moa_E_72
#> 20: G00005_drug_005_moa_A_G00026_drug_026_moa_E_72
#> 21: G00006_drug_006_moa_A_G00021_drug_021_moa_D_72
#> 22: G00006_drug_006_moa_A_G00021_drug_021_moa_D_72
#> 23: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 24: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#>                                               r Id
#>                                             <char>
#>                                c Id Normalization Type Fit Source Bliss Score
#>                              <char>             <char>     <char>       <num>
#>  1: CL00016_cellline_GB_tissue_y_46                 GR        gDR 0.039386959
#>  2: CL00016_cellline_GB_tissue_y_46                 RV        gDR 0.011889039
#>  3: CL00016_cellline_GB_tissue_y_46                 GR        gDR 0.045591354
#>  4: CL00016_cellline_GB_tissue_y_46                 RV        gDR 0.013056661
#>  5: CL00016_cellline_GB_tissue_y_46                 GR        gDR 0.045176579
#>  6: CL00016_cellline_GB_tissue_y_46                 RV        gDR 0.014644613
#>  7: CL00016_cellline_GB_tissue_y_46                 GR        gDR 0.046993939
#>  8: CL00016_cellline_GB_tissue_y_46                 RV        gDR 0.016299122
#>  9: CL00016_cellline_GB_tissue_y_46                 GR        gDR 0.020293690
#> 10: CL00016_cellline_GB_tissue_y_46                 RV        gDR 0.005337033
#> 11: CL00016_cellline_GB_tissue_y_46                 GR        gDR 0.022952967
#> 12: CL00016_cellline_GB_tissue_y_46                 RV        gDR 0.005529326
#> 13: CL00017_cellline_HB_tissue_y_50                 GR        gDR 0.024430399
#> 14: CL00017_cellline_HB_tissue_y_50                 RV        gDR 0.008970004
#> 15: CL00017_cellline_HB_tissue_y_50                 GR        gDR 0.047794901
#> 16: CL00017_cellline_HB_tissue_y_50                 RV        gDR 0.014493716
#> 17: CL00017_cellline_HB_tissue_y_50                 GR        gDR 0.032388931
#> 18: CL00017_cellline_HB_tissue_y_50                 RV        gDR 0.013518264
#> 19: CL00017_cellline_HB_tissue_y_50                 GR        gDR 0.048273626
#> 20: CL00017_cellline_HB_tissue_y_50                 RV        gDR 0.014159232
#> 21: CL00017_cellline_HB_tissue_y_50                 GR        gDR 0.077379460
#> 22: CL00017_cellline_HB_tissue_y_50                 RV        gDR 0.029531491
#> 23: CL00017_cellline_HB_tissue_y_50                 GR        gDR 0.139127666
#> 24: CL00017_cellline_HB_tissue_y_50                 RV        gDR 0.056935611
#>                                c Id Normalization Type Fit Source Bliss Score
#>                              <char>             <char>     <char>       <num>
#>       HSA Score  CIScore 50   CIScore 80 Gnumber Drug Name Drug MOA Gnumber 2
#>           <num>       <num>        <num>  <char>    <char>   <char>    <char>
#>  1: 0.039386959  0.01747735  0.029337758  G00004  drug_004    moa_A    G00021
#>  2: 0.011889039  0.03485816  0.042184690  G00004  drug_004    moa_A    G00021
#>  3: 0.045591354  0.02293000  0.010300148  G00004  drug_004    moa_A    G00026
#>  4: 0.013056661  0.01554275  0.013571845  G00004  drug_004    moa_A    G00026
#>  5: 0.045176579  0.01476898  0.027162472  G00005  drug_005    moa_A    G00021
#>  6: 0.014644613  0.02046719  0.032236614  G00005  drug_005    moa_A    G00021
#>  7: 0.046993939  0.15986726 -0.048280684  G00005  drug_005    moa_A    G00026
#>  8: 0.016299122  0.01524069  0.022185871  G00005  drug_005    moa_A    G00026
#>  9: 0.020293690  0.02309152  0.027147952  G00006  drug_006    moa_A    G00021
#> 10: 0.005337033  0.02563107  0.021081550  G00006  drug_006    moa_A    G00021
#> 11: 0.022952967  0.02039404  0.019774471  G00006  drug_006    moa_A    G00026
#> 12: 0.005529326  0.02106643  0.020332834  G00006  drug_006    moa_A    G00026
#> 13: 0.024430399  0.02002235  0.015382529  G00004  drug_004    moa_A    G00021
#> 14: 0.008970004  0.01573286  0.014713425  G00004  drug_004    moa_A    G00021
#> 15: 0.047794901  0.02639060  0.021365112  G00004  drug_004    moa_A    G00026
#> 16: 0.014493716  0.02122847  0.016280462  G00004  drug_004    moa_A    G00026
#> 17: 0.032388931  0.02205106 -0.005903405  G00005  drug_005    moa_A    G00021
#> 18: 0.013518264 -0.00401062 -0.032211087  G00005  drug_005    moa_A    G00021
#> 19: 0.048273626  0.03074943  0.051788631  G00005  drug_005    moa_A    G00026
#> 20: 0.014159232  0.04725145  0.081596255  G00005  drug_005    moa_A    G00026
#> 21: 0.077379460 -0.01169389 -0.138659682  G00006  drug_006    moa_A    G00021
#> 22: 0.029531491 -0.07492990 -0.297592786  G00006  drug_006    moa_A    G00021
#> 23: 0.139127666  0.07855021  0.339741093  G00006  drug_006    moa_A    G00026
#> 24: 0.056935611  0.16251115  0.384075521  G00006  drug_006    moa_A    G00026
#>       HSA Score  CIScore 50   CIScore 80 Gnumber Drug Name Drug MOA Gnumber 2
#>           <num>       <num>        <num>  <char>    <char>   <char>    <char>
#>     Drug Name 2 Drug MOA 2 Duration    Clid Cell Line Name   Tissue
#>          <char>     <char>    <num>  <char>         <char>   <char>
#>  1:    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#>  2:    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#>  3:    drug_026      moa_E       72 CL00016    cellline_GB tissue_y
#>  4:    drug_026      moa_E       72 CL00016    cellline_GB tissue_y
#>  5:    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#>  6:    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#>  7:    drug_026      moa_E       72 CL00016    cellline_GB tissue_y
#>  8:    drug_026      moa_E       72 CL00016    cellline_GB tissue_y
#>  9:    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#> 10:    drug_021      moa_D       72 CL00016    cellline_GB tissue_y
#> 11:    drug_026      moa_E       72 CL00016    cellline_GB tissue_y
#> 12:    drug_026      moa_E       72 CL00016    cellline_GB tissue_y
#> 13:    drug_021      moa_D       72 CL00017    cellline_HB tissue_y
#> 14:    drug_021      moa_D       72 CL00017    cellline_HB tissue_y
#> 15:    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#> 16:    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#> 17:    drug_021      moa_D       72 CL00017    cellline_HB tissue_y
#> 18:    drug_021      moa_D       72 CL00017    cellline_HB tissue_y
#> 19:    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#> 20:    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#> 21:    drug_021      moa_D       72 CL00017    cellline_HB tissue_y
#> 22:    drug_021      moa_D       72 CL00017    cellline_HB tissue_y
#> 23:    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#> 24:    drug_026      moa_E       72 CL00017    cellline_HB tissue_y
#>     Drug Name 2 Drug MOA 2 Duration    Clid Cell Line Name   Tissue
#>          <char>     <char>    <num>  <char>         <char>   <char>
#>     Reference Division Time
#>                       <num>
#>  1:                      46
#>  2:                      46
#>  3:                      46
#>  4:                      46
#>  5:                      46
#>  6:                      46
#>  7:                      46
#>  8:                      46
#>  9:                      46
#> 10:                      46
#> 11:                      46
#> 12:                      46
#> 13:                      50
#> 14:                      50
#> 15:                      50
#> 16:                      50
#> 17:                      50
#> 18:                      50
#> 19:                      50
#> 20:                      50
#> 21:                      50
#> 22:                      50
#> 23:                      50
#> 24:                      50
#>     Reference Division Time
#>                       <num>
#> 
#> $isobolograms
#>                                                  r Id
#>                                                <char>
#>     1: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>     2: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>     3: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>     4: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>     5: G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#>    ---                                               
#> 12948: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 12949: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 12950: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 12951: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 12952: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#>                                   c Id       Pos x     Pos y   Pos x Ref
#>                                 <char>       <num>     <num>       <num>
#>     1: CL00016_cellline_GB_tissue_y_46 -3.67577834 -2.174830 -3.67577834
#>     2: CL00016_cellline_GB_tissue_y_46 -3.16068978 -2.170894 -3.16465561
#>     3: CL00016_cellline_GB_tissue_y_46 -3.01855302 -2.170179 -3.02324569
#>     4: CL00016_cellline_GB_tissue_y_46 -2.87575556 -2.168803 -2.88184018
#>     5: CL00016_cellline_GB_tissue_y_46 -2.73307544 -2.167544 -2.74044077
#>    ---                                                                  
#> 12948: CL00017_cellline_HB_tissue_y_50  0.07749422 -2.068615  0.03788669
#> 12949: CL00017_cellline_HB_tissue_y_50  0.73824771 -2.843478  0.73824771
#> 12950: CL00017_cellline_HB_tissue_y_50  0.73824771 -2.984900  0.73824771
#> 12951: CL00017_cellline_HB_tissue_y_50  0.73824771 -3.126321  0.73824771
#> 12952: CL00017_cellline_HB_tissue_y_50  0.73824771 -3.675778  0.73824771
#>        Pos y Ref Log10 Ratio Conc    Log2 CI Iso Level     CI 100x
#>            <num>            <num>      <num>    <char>       <num>
#>     1: -2.174830       -1.0613306 0.00000000     -0.25  0.03872664
#>     2: -2.174860       -0.6998913 0.01224389     -0.25  0.03872664
#>     3: -2.174871       -0.5998913 0.01836066     -0.25  0.03872664
#>     4: -2.174887       -0.4998913 0.02281227     -0.25  0.03872664
#>     5: -2.174909       -0.3998913 0.02902613     -0.25  0.03872664
#>    ---                                                            
#> 12948: -2.108222        1.5175283 0.14930479       0.8 -0.02939538
#> 12949: -2.909123        2.7175283 0.15685335       0.8 -0.02939538
#> 12950: -3.033195        2.8175283 0.12396001       0.8 -0.02939538
#> 12951: -3.161643        2.9175283 0.00000000       0.8 -0.02939538
#> 12952: -3.675778        3.3062746 0.00000000       0.8 -0.02939538
#>        Normalization Type Gnumber Drug Name Drug MOA Gnumber 2 Drug Name 2
#>                    <char>  <char>    <char>   <char>    <char>      <char>
#>     1:                 GR  G00004  drug_004    moa_A    G00021    drug_021
#>     2:                 GR  G00004  drug_004    moa_A    G00021    drug_021
#>     3:                 GR  G00004  drug_004    moa_A    G00021    drug_021
#>     4:                 GR  G00004  drug_004    moa_A    G00021    drug_021
#>     5:                 GR  G00004  drug_004    moa_A    G00021    drug_021
#>    ---                                                                    
#> 12948:                 RV  G00006  drug_006    moa_A    G00026    drug_026
#> 12949:                 RV  G00006  drug_006    moa_A    G00026    drug_026
#> 12950:                 RV  G00006  drug_006    moa_A    G00026    drug_026
#> 12951:                 RV  G00006  drug_006    moa_A    G00026    drug_026
#> 12952:                 RV  G00006  drug_006    moa_A    G00026    drug_026
#>        Drug MOA 2 Duration    Clid Cell Line Name   Tissue
#>            <char>    <num>  <char>         <char>   <char>
#>     1:      moa_D       72 CL00016    cellline_GB tissue_y
#>     2:      moa_D       72 CL00016    cellline_GB tissue_y
#>     3:      moa_D       72 CL00016    cellline_GB tissue_y
#>     4:      moa_D       72 CL00016    cellline_GB tissue_y
#>     5:      moa_D       72 CL00016    cellline_GB tissue_y
#>    ---                                                    
#> 12948:      moa_E       72 CL00017    cellline_HB tissue_y
#> 12949:      moa_E       72 CL00017    cellline_HB tissue_y
#> 12950:      moa_E       72 CL00017    cellline_HB tissue_y
#> 12951:      moa_E       72 CL00017    cellline_HB tissue_y
#> 12952:      moa_E       72 CL00017    cellline_HB tissue_y
#>        Reference Division Time
#>                          <num>
#>     1:                      46
#>     2:                      46
#>     3:                      46
#>     4:                      46
#>     5:                      46
#>    ---                        
#> 12948:                      50
#> 12949:                      50
#> 12950:                      50
#> 12951:                      50
#> 12952:                      50
#> 
```
