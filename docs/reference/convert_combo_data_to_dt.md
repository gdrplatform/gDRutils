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
#>                                  c Id Concentration Concentration 2     Smooth
#>                                <char>         <num>           <num>      <num>
#>    1: CL00016_cellline_GB_tissue_y_46          0.00         0.00000  1.0000000
#>    2: CL00016_cellline_GB_tissue_y_46          0.00         0.00000  1.0000000
#>    3: CL00016_cellline_GB_tissue_y_46          0.00         0.00100  1.0000000
#>    4: CL00016_cellline_GB_tissue_y_46          0.00         0.00100  1.0000000
#>    5: CL00016_cellline_GB_tissue_y_46          0.00         0.00316  1.0000000
#>   ---                                                                         
#> 1940: CL00017_cellline_HB_tissue_y_50          3.16         0.31600  0.1828034
#> 1941: CL00017_cellline_HB_tissue_y_50          3.16         1.00000 -0.4035902
#> 1942: CL00017_cellline_HB_tissue_y_50          3.16         1.00000  0.1918578
#> 1943: CL00017_cellline_HB_tissue_y_50          3.16         3.16000 -0.4008597
#> 1944: CL00017_cellline_HB_tissue_y_50          3.16         3.16000  0.1921170
#>       Normalization Type   HSA Excess Bliss Excess Gnumber Drug Name Drug MOA
#>                   <char>        <num>        <num>  <char>    <char>   <char>
#>    1:                 GR           NA           NA  G00004  drug_004    moa_A
#>    2:                 RV           NA           NA  G00004  drug_004    moa_A
#>    3:                 GR           NA           NA  G00004  drug_004    moa_A
#>    4:                 RV           NA           NA  G00004  drug_004    moa_A
#>    5:                 GR           NA           NA  G00004  drug_004    moa_A
#>   ---                                                                        
#> 1940:                 RV -0.001579471 -0.001579471  G00006  drug_006    moa_A
#> 1941:                 GR  0.003961392  0.003961392  G00006  drug_006    moa_A
#> 1942:                 RV -0.010633816 -0.010633816  G00006  drug_006    moa_A
#> 1943:                 GR  0.001230976  0.001230976  G00006  drug_006    moa_A
#> 1944:                 RV -0.010893026 -0.010893026  G00006  drug_006    moa_A
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
#>                                c Id Normalization Type    HSA Score
#>                              <char>             <char>        <num>
#>  1: CL00016_cellline_GB_tissue_y_46                 GR 0.0191864677
#>  2: CL00016_cellline_GB_tissue_y_46                 RV 0.0054825946
#>  3: CL00016_cellline_GB_tissue_y_46                 GR 0.0008640781
#>  4: CL00016_cellline_GB_tissue_y_46                 RV 0.0002808689
#>  5: CL00016_cellline_GB_tissue_y_46                 GR 0.0212207909
#>  6: CL00016_cellline_GB_tissue_y_46                 RV 0.0064007236
#>  7: CL00016_cellline_GB_tissue_y_46                 GR 0.0282311962
#>  8: CL00016_cellline_GB_tissue_y_46                 RV 0.0004154646
#>  9: CL00016_cellline_GB_tissue_y_46                 GR 0.0142686219
#> 10: CL00016_cellline_GB_tissue_y_46                 RV 0.0035563938
#> 11: CL00016_cellline_GB_tissue_y_46                 GR 0.0006824547
#> 12: CL00016_cellline_GB_tissue_y_46                 RV 0.0002571149
#> 13: CL00017_cellline_HB_tissue_y_50                 GR 0.0195462835
#> 14: CL00017_cellline_HB_tissue_y_50                 RV 0.0062733020
#> 15: CL00017_cellline_HB_tissue_y_50                 GR 0.0047716696
#> 16: CL00017_cellline_HB_tissue_y_50                 RV 0.0020015238
#> 17: CL00017_cellline_HB_tissue_y_50                 GR 0.0241474816
#> 18: CL00017_cellline_HB_tissue_y_50                 RV 0.0080924532
#> 19: CL00017_cellline_HB_tissue_y_50                 GR 0.0084478065
#> 20: CL00017_cellline_HB_tissue_y_50                 RV 0.0009734877
#> 21: CL00017_cellline_HB_tissue_y_50                 GR 0.0550440010
#> 22: CL00017_cellline_HB_tissue_y_50                 RV 0.0268261867
#> 23: CL00017_cellline_HB_tissue_y_50                 GR 0.0344185361
#> 24: CL00017_cellline_HB_tissue_y_50                 RV 0.0090907653
#>                                c Id Normalization Type    HSA Score
#>                              <char>             <char>        <num>
#>      Bliss Score  CIScore 50   CIScore 80 Gnumber Drug Name Drug MOA Gnumber 2
#>            <num>       <num>        <num>  <char>    <char>   <char>    <char>
#>  1: 0.0191864677  0.01747735  0.029337758  G00004  drug_004    moa_A    G00021
#>  2: 0.0054825946  0.03485816  0.042184690  G00004  drug_004    moa_A    G00021
#>  3: 0.0008640781  0.02292997  0.010300132  G00004  drug_004    moa_A    G00026
#>  4: 0.0002808689  0.01554275  0.013571845  G00004  drug_004    moa_A    G00026
#>  5: 0.0212207909  0.01531961  0.027161987  G00005  drug_005    moa_A    G00021
#>  6: 0.0064007236  0.02046719  0.032236614  G00005  drug_005    moa_A    G00021
#>  7: 0.0282311962  0.16264958 -0.048281369  G00005  drug_005    moa_A    G00026
#>  8: 0.0004154646  0.01523954  0.022184526  G00005  drug_005    moa_A    G00026
#>  9: 0.0142686219  0.02309152  0.027147952  G00006  drug_006    moa_A    G00021
#> 10: 0.0035563938  0.02563107  0.021081550  G00006  drug_006    moa_A    G00021
#> 11: 0.0006824547  0.02039404  0.019774471  G00006  drug_006    moa_A    G00026
#> 12: 0.0002571149  0.02106643  0.020332834  G00006  drug_006    moa_A    G00026
#> 13: 0.0195462835  0.02002235  0.015382530  G00004  drug_004    moa_A    G00021
#> 14: 0.0062733020  0.01573286  0.014713425  G00004  drug_004    moa_A    G00021
#> 15: 0.0047716696  0.02639060  0.021365112  G00004  drug_004    moa_A    G00026
#> 16: 0.0020015238  0.02122847  0.016280462  G00004  drug_004    moa_A    G00026
#> 17: 0.0241474816  0.02205106 -0.005903405  G00005  drug_005    moa_A    G00021
#> 18: 0.0080924532 -0.00401062 -0.032211087  G00005  drug_005    moa_A    G00021
#> 19: 0.0084478065  0.03074943  0.051788631  G00005  drug_005    moa_A    G00026
#> 20: 0.0009734877  0.04725145  0.081596255  G00005  drug_005    moa_A    G00026
#> 21: 0.0550440010 -0.01169389 -0.138659682  G00006  drug_006    moa_A    G00021
#> 22: 0.0268261867 -0.07492990 -0.297592786  G00006  drug_006    moa_A    G00021
#> 23: 0.0344185361  0.07855021  0.339741093  G00006  drug_006    moa_A    G00026
#> 24: 0.0090907653  0.16251115  0.384075521  G00006  drug_006    moa_A    G00026
#>      Bliss Score  CIScore 50   CIScore 80 Gnumber Drug Name Drug MOA Gnumber 2
#>            <num>       <num>        <num>  <char>    <char>   <char>    <char>
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
#> 12990: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 12991: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 12992: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 12993: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#> 12994: G00006_drug_006_moa_A_G00026_drug_026_moa_E_72
#>                                   c Id Iso Level      Pos x     Pos y
#>                                 <char>    <char>      <num>     <num>
#>     1: CL00016_cellline_GB_tissue_y_46     -0.05 -3.6757783 -2.283728
#>     2: CL00016_cellline_GB_tissue_y_46     -0.05 -3.1612132 -2.280315
#>     3: CL00016_cellline_GB_tissue_y_46     -0.05 -3.0184836 -2.279007
#>     4: CL00016_cellline_GB_tissue_y_46     -0.05 -2.8761652 -2.278110
#>     5: CL00016_cellline_GB_tissue_y_46     -0.05 -2.7335957 -2.276962
#>    ---                                                               
#> 12990: CL00017_cellline_HB_tissue_y_50       0.8  0.7382477 -3.126321
#> 12991: CL00017_cellline_HB_tissue_y_50       0.8  0.7382477 -3.000234
#> 12992: CL00017_cellline_HB_tissue_y_50       0.8  0.7382477 -2.984900
#> 12993: CL00017_cellline_HB_tissue_y_50       0.8  0.7382477 -2.858813
#> 12994: CL00017_cellline_HB_tissue_y_50       0.8  0.7382477 -2.843478
#>         Pos x Ref Pos y Ref Log10 Ratio Conc    Log2 CI Normalization Type
#>             <num>     <num>            <num>      <num>             <char>
#>     1: -3.6757783 -2.283728       -0.9843283 0.00000000                 GR
#>     2: -3.1646556 -2.283758       -0.6228889 0.01152660                 GR
#>     3: -3.0232457 -2.283769       -0.5228889 0.01721126                 GR
#>     4: -2.8818402 -2.283785       -0.4228889 0.02231182                 GR
#>     5: -2.7404408 -2.283807       -0.3228889 0.03271470                 GR
#>    ---                                                                    
#> 12990:  0.7382477 -3.161643        2.9175283 0.00000000                 RV
#> 12991:  0.7382477 -3.065879        2.8283715 0.16819143                 GR
#> 12992:  0.7382477 -3.033195        2.8175283 0.12396001                 RV
#> 12993:  0.7382477 -2.947395        2.7283715 0.14361810                 GR
#> 12994:  0.7382477 -2.909123        2.7175283 0.15685335                 RV
#>        Gnumber Drug Name Drug MOA Gnumber 2 Drug Name 2 Drug MOA 2 Duration
#>         <char>    <char>   <char>    <char>      <char>     <char>    <num>
#>     1:  G00004  drug_004    moa_A    G00021    drug_021      moa_D       72
#>     2:  G00004  drug_004    moa_A    G00021    drug_021      moa_D       72
#>     3:  G00004  drug_004    moa_A    G00021    drug_021      moa_D       72
#>     4:  G00004  drug_004    moa_A    G00021    drug_021      moa_D       72
#>     5:  G00004  drug_004    moa_A    G00021    drug_021      moa_D       72
#>    ---                                                                     
#> 12990:  G00006  drug_006    moa_A    G00026    drug_026      moa_E       72
#> 12991:  G00006  drug_006    moa_A    G00026    drug_026      moa_E       72
#> 12992:  G00006  drug_006    moa_A    G00026    drug_026      moa_E       72
#> 12993:  G00006  drug_006    moa_A    G00026    drug_026      moa_E       72
#> 12994:  G00006  drug_006    moa_A    G00026    drug_026      moa_E       72
#>           Clid Cell Line Name   Tissue Reference Division Time
#>         <char>         <char>   <char>                   <num>
#>     1: CL00016    cellline_GB tissue_y                      46
#>     2: CL00016    cellline_GB tissue_y                      46
#>     3: CL00016    cellline_GB tissue_y                      46
#>     4: CL00016    cellline_GB tissue_y                      46
#>     5: CL00016    cellline_GB tissue_y                      46
#>    ---                                                        
#> 12990: CL00017    cellline_HB tissue_y                      50
#> 12991: CL00017    cellline_HB tissue_y                      50
#> 12992: CL00017    cellline_HB tissue_y                      50
#> 12993: CL00017    cellline_HB tissue_y                      50
#> 12994: CL00017    cellline_HB tissue_y                      50
#> 
```
