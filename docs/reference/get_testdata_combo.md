# get_testdata_combo

Function to obtain data from gDRtestData and prepare for unit tests

## Usage

``` r
get_testdata_combo()
```

## Value

list with drugs, cell_lines, raw_data and assay_data

## Examples

``` r
get_testdata_combo()
#> $drug_names
#> [1] "drug_001" "drug_002" "drug_011" "drug_021" "drug_026" "drug_031"
#> 
#> $cell_line_names
#> [1] "cellline_AA" "cellline_EA" "cellline_IB" "cellline_MC" "cellline_BC"
#> [6] "cellline_FD" "cellline_JE" "cellline_NE"
#> 
#> $dt
#>     Drug Name Cell Line Name                                            rId
#>        <char>         <char>                                         <char>
#>  1:  drug_001    cellline_AA G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#>  2:  drug_002    cellline_AA G00002_drug_002_moa_A_G00021_drug_021_moa_D_72
#>  3:  drug_011    cellline_AA G00011_drug_011_moa_B_G00021_drug_021_moa_D_72
#>  4:  drug_001    cellline_EA G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#>  5:  drug_002    cellline_EA G00002_drug_002_moa_A_G00021_drug_021_moa_D_72
#>  6:  drug_011    cellline_EA G00011_drug_011_moa_B_G00021_drug_021_moa_D_72
#>  7:  drug_001    cellline_IB G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#>  8:  drug_002    cellline_IB G00002_drug_002_moa_A_G00021_drug_021_moa_D_72
#>  9:  drug_011    cellline_IB G00011_drug_011_moa_B_G00021_drug_021_moa_D_72
#> 10:  drug_001    cellline_MC G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#> 11:  drug_002    cellline_MC G00002_drug_002_moa_A_G00021_drug_021_moa_D_72
#> 12:  drug_011    cellline_MC G00011_drug_011_moa_B_G00021_drug_021_moa_D_72
#> 13:  drug_001    cellline_BC G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#> 14:  drug_002    cellline_BC G00002_drug_002_moa_A_G00021_drug_021_moa_D_72
#> 15:  drug_011    cellline_BC G00011_drug_011_moa_B_G00021_drug_021_moa_D_72
#> 16:  drug_001    cellline_FD G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#> 17:  drug_002    cellline_FD G00002_drug_002_moa_A_G00021_drug_021_moa_D_72
#> 18:  drug_011    cellline_FD G00011_drug_011_moa_B_G00021_drug_021_moa_D_72
#> 19:  drug_001    cellline_JE G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#> 20:  drug_002    cellline_JE G00002_drug_002_moa_A_G00021_drug_021_moa_D_72
#> 21:  drug_011    cellline_JE G00011_drug_011_moa_B_G00021_drug_021_moa_D_72
#> 22:  drug_001    cellline_NE G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#> 23:  drug_002    cellline_NE G00002_drug_002_moa_A_G00021_drug_021_moa_D_72
#> 24:  drug_011    cellline_NE G00011_drug_011_moa_B_G00021_drug_021_moa_D_72
#> 25:  drug_021    cellline_AA                       G00021_drug_021_moa_D_72
#> 26:  drug_026    cellline_AA                       G00026_drug_026_moa_E_72
#> 27:  drug_031    cellline_AA                       G00031_drug_031_moa_F_72
#> 28:  drug_021    cellline_EA                       G00021_drug_021_moa_D_72
#> 29:  drug_026    cellline_EA                       G00026_drug_026_moa_E_72
#> 30:  drug_031    cellline_EA                       G00031_drug_031_moa_F_72
#> 31:  drug_021    cellline_IB                       G00021_drug_021_moa_D_72
#> 32:  drug_026    cellline_IB                       G00026_drug_026_moa_E_72
#> 33:  drug_031    cellline_IB                       G00031_drug_031_moa_F_72
#> 34:  drug_021    cellline_MC                       G00021_drug_021_moa_D_72
#> 35:  drug_026    cellline_MC                       G00026_drug_026_moa_E_72
#> 36:  drug_031    cellline_MC                       G00031_drug_031_moa_F_72
#> 37:  drug_021    cellline_BC                       G00021_drug_021_moa_D_72
#> 38:  drug_026    cellline_BC                       G00026_drug_026_moa_E_72
#> 39:  drug_031    cellline_BC                       G00031_drug_031_moa_F_72
#> 40:  drug_021    cellline_FD                       G00021_drug_021_moa_D_72
#> 41:  drug_026    cellline_FD                       G00026_drug_026_moa_E_72
#> 42:  drug_031    cellline_FD                       G00031_drug_031_moa_F_72
#> 43:  drug_021    cellline_JE                       G00021_drug_021_moa_D_72
#> 44:  drug_026    cellline_JE                       G00026_drug_026_moa_E_72
#> 45:  drug_031    cellline_JE                       G00031_drug_031_moa_F_72
#> 46:  drug_021    cellline_NE                       G00021_drug_021_moa_D_72
#> 47:  drug_026    cellline_NE                       G00026_drug_026_moa_E_72
#> 48:  drug_031    cellline_NE                       G00031_drug_031_moa_F_72
#>     Drug Name Cell Line Name                                            rId
#>        <char>         <char>                                         <char>
#>                                 cId        EC50         x_AOC          h RV
#>                              <char>       <num>         <num>         <num>
#>  1: CL00010_cellline_AA_tissue_x_22  0.31359259  0.6864074074  0.6864074074
#>  2: CL00010_cellline_AA_tissue_x_22  0.46488148  0.5351185185  0.5351185185
#>  3: CL00010_cellline_AA_tissue_x_22  0.30005185  0.6999481481  0.6999481481
#>  4: CL00014_cellline_EA_tissue_x_38  0.54018519  0.4598148148  0.4598148148
#>  5: CL00014_cellline_EA_tissue_x_38  0.46981481  0.5301851852  0.5301851852
#>  6: CL00014_cellline_EA_tissue_x_38 -0.27716296  1.2771629630  1.2771629630
#>  7: CL00018_cellline_IB_tissue_y_54 -0.70610370  1.7061037037  1.7061037037
#>  8: CL00018_cellline_IB_tissue_y_54 -0.86958519  1.8695851852  1.8695851852
#>  9: CL00018_cellline_IB_tissue_y_54 -0.60236296  1.6023629630  1.6023629630
#> 10: CL00022_cellline_MC_tissue_z_70 -0.73374444  1.7337444444  1.7337444444
#> 11: CL00022_cellline_MC_tissue_z_70 -0.58530000  1.5853000000  1.5853000000
#> 12: CL00022_cellline_MC_tissue_z_70 -0.54707037  1.5470703704  1.5470703704
#> 13: CL00026_cellline_BC_tissue_w_26 -0.27410000  1.2741000000  1.2741000000
#> 14: CL00026_cellline_BC_tissue_w_26  0.19546296  0.8045370370  0.8045370370
#> 15: CL00026_cellline_BC_tissue_w_26  0.08634815  0.9136518519  0.9136518519
#> 16: CL00030_cellline_FD_tissue_w_42 -0.36797407  1.3679740741  1.3679740741
#> 17: CL00030_cellline_FD_tissue_w_42 -0.30042222  1.3004222222  1.3004222222
#> 18: CL00030_cellline_FD_tissue_w_42 -0.41170741  1.4117074074  1.4117074074
#> 19: CL00034_cellline_JE_tissue_w_58 -0.69018148  1.6901814815  1.6901814815
#> 20: CL00034_cellline_JE_tissue_w_58 -0.64542222  1.6454222222  1.6454222222
#> 21: CL00034_cellline_JE_tissue_w_58 -0.64786296  1.6478629630  1.6478629630
#> 22: CL00038_cellline_NE_tissue_w_74 -0.75242963  1.7524296296  1.7524296296
#> 23: CL00038_cellline_NE_tissue_w_74 -0.54112963  1.5411296296  1.5411296296
#> 24: CL00038_cellline_NE_tissue_w_74 -0.43574815  1.4357481481  1.4357481481
#> 25: CL00010_cellline_AA_tissue_x_22  0.99633210  0.0036679012  0.0036679012
#> 26: CL00010_cellline_AA_tissue_x_22  0.99605432  0.0039456790  0.0039456790
#> 27: CL00010_cellline_AA_tissue_x_22  0.99742531  0.0025746914  0.0025746914
#> 28: CL00014_cellline_EA_tissue_x_38  0.99882160  0.0011783951  0.0011783951
#> 29: CL00014_cellline_EA_tissue_x_38  1.00082469 -0.0008246914 -0.0008246914
#> 30: CL00014_cellline_EA_tissue_x_38  0.99524444  0.0047555556  0.0047555556
#> 31: CL00018_cellline_IB_tissue_y_54  0.99930432  0.0006956790  0.0006956790
#> 32: CL00018_cellline_IB_tissue_y_54  0.99741975  0.0025802469  0.0025802469
#> 33: CL00018_cellline_IB_tissue_y_54  0.99455556  0.0054444444  0.0054444444
#> 34: CL00022_cellline_MC_tissue_z_70  0.99650309  0.0034969136  0.0034969136
#> 35: CL00022_cellline_MC_tissue_z_70  1.00027901 -0.0002790123 -0.0002790123
#> 36: CL00022_cellline_MC_tissue_z_70  0.99801481  0.0019851852  0.0019851852
#> 37: CL00026_cellline_BC_tissue_w_26  0.99741235  0.0025876543  0.0025876543
#> 38: CL00026_cellline_BC_tissue_w_26  1.00112284 -0.0011228395 -0.0011228395
#> 39: CL00026_cellline_BC_tissue_w_26  0.99772222  0.0022777778  0.0022777778
#> 40: CL00030_cellline_FD_tissue_w_42  0.99697037  0.0030296296  0.0030296296
#> 41: CL00030_cellline_FD_tissue_w_42  1.00156420 -0.0015641975 -0.0015641975
#> 42: CL00030_cellline_FD_tissue_w_42  0.99700185  0.0029981481  0.0029981481
#> 43: CL00034_cellline_JE_tissue_w_58  1.00209074 -0.0020907407 -0.0020907407
#> 44: CL00034_cellline_JE_tissue_w_58  1.00434630 -0.0043462963 -0.0043462963
#> 45: CL00034_cellline_JE_tissue_w_58  1.00456852 -0.0045685185 -0.0045685185
#> 46: CL00038_cellline_NE_tissue_w_74  0.99237531  0.0076246914  0.0076246914
#> 47: CL00038_cellline_NE_tissue_w_74  0.99570370  0.0042962963  0.0042962963
#> 48: CL00038_cellline_NE_tissue_w_74  0.99102160  0.0089783951  0.0089783951
#>                                 cId        EC50         x_AOC          h RV
#>                              <char>       <num>         <num>         <num>
#>     GEC50      GR Max  ec50      GR Inf        GR 0  h GR E Inf   p_value
#>     <num>       <num> <num>       <num>       <num> <num> <num>     <num>
#>  1:  -Inf  0.31106667     0  0.31359259  0.31359259 1e-04     0 1.0000000
#>  2:  -Inf  0.44736667     0  0.46488148  0.46488148 1e-04     0 0.8071221
#>  3:  -Inf  0.28696667     0  0.30005185  0.30005185 1e-04     0 0.8857690
#>  4:   Inf  0.53176667     0  0.54018519  0.54018519 1e-04     0 0.9931766
#>  5:  -Inf  0.45896667     0  0.46981481  0.46981481 1e-04     0 1.0000000
#>  6:  -Inf -0.31976667     0 -0.27716296 -0.27716296 1e-04     0 0.7207726
#>  7:  -Inf -0.72603333     0 -0.70610370 -0.70610370 1e-04     0 0.9103471
#>  8:  -Inf -0.91793333     0 -0.86958519 -0.86958519 1e-04     0 1.0000000
#>  9:  -Inf -0.61163333     0 -0.60236296 -0.60236296 1e-04     0 1.0000000
#> 10:  -Inf -0.74943333     0 -0.73374444 -0.73374444 1e-04     0 0.5342941
#> 11:  -Inf -0.60303333     0 -0.58530000 -0.58530000 1e-04     0 1.0000000
#> 12:  -Inf -0.52776667     0 -0.54707037 -0.54707037 1e-04     0 1.0000000
#> 13:  -Inf -0.26623333     0 -0.27410000 -0.27410000 1e-04     0 0.9488407
#> 14:  -Inf  0.19010000     0  0.19546296  0.19546296 1e-04     0 1.0000000
#> 15:  -Inf  0.08556667     0  0.08634815  0.08634815 1e-04     0 0.7266117
#> 16:  -Inf -0.36960000     0 -0.36797407 -0.36797407 1e-04     0 1.0000000
#> 17:  -Inf -0.28823333     0 -0.30042222 -0.30042222 1e-04     0 0.9333698
#> 18:  -Inf -0.42910000     0 -0.41170741 -0.41170741 1e-04     0 0.9975568
#> 19:  -Inf -0.70030000     0 -0.69018148 -0.69018148 1e-04     0 1.0000000
#> 20:  -Inf -0.65620000     0 -0.64542222 -0.64542222 1e-04     0 0.2453094
#> 21:  -Inf -0.64646667     0 -0.64786296 -0.64786296 1e-04     0 1.0000000
#> 22:  -Inf -0.77066667     0 -0.75242963 -0.75242963 1e-04     0 0.8325653
#> 23:  -Inf -0.54173333     0 -0.54112963 -0.54112963 1e-04     0 1.0000000
#> 24:  -Inf -0.43120000     0 -0.43574815 -0.43574815 1e-04     0 0.3266012
#> 25:   Inf  0.98825556     0  0.99633210  0.99633210 1e-04     0 0.9083469
#> 26:   Inf  0.98951667     0  0.99605432  0.99605432 1e-04     0 0.5747625
#> 27:   Inf  0.99997222     0  0.99742531  0.99742531 1e-04     0 1.0000000
#> 28:   Inf  1.00171111     0  0.99882160  0.99882160 1e-04     0 1.0000000
#> 29:   Inf  0.99019444     0  1.00082469  1.00082469 1e-04     0 1.0000000
#> 30:   Inf  0.98635000     0  0.99524444  0.99524444 1e-04     0 0.7599590
#> 31:   Inf  0.99443333     0  0.99930432  0.99930432 1e-04     0 1.0000000
#> 32:   Inf  0.98379444     0  0.99741975  0.99741975 1e-04     0 0.8885349
#> 33:   Inf  0.99092222     0  0.99455556  0.99455556 1e-04     0 0.8825499
#> 34:   Inf  0.99152778     0  0.99650309  0.99650309 1e-04     0 0.8933194
#> 35:   Inf  1.00274444     0  1.00027901  1.00027901 1e-04     0 1.0000000
#> 36:   Inf  0.99694444     0  0.99801481  0.99801481 1e-04     0 1.0000000
#> 37:   Inf  0.99417222     0  0.99741235  0.99741235 1e-04     0 1.0000000
#> 38:   Inf  0.99040000     0  1.00112284  1.00112284 1e-04     0 1.0000000
#> 39:   Inf  0.98581667     0  0.99772222  0.99772222 1e-04     0 1.0000000
#> 40:   Inf  0.99617778     0  0.99697037  0.99697037 1e-04     0 0.9235286
#> 41:   Inf  0.99464444     0  1.00156420  1.00156420 1e-04     0 1.0000000
#> 42:   Inf  0.99352222     0  0.99700185  0.99700185 1e-04     0 1.0000000
#> 43:   Inf  0.99495556     0  1.00209074  1.00209074 1e-04     0 1.0000000
#> 44:   Inf  0.99893889     0  1.00434630  1.00434630 1e-04     0 1.0000000
#> 45:   Inf  0.99412222     0  1.00456852  1.00456852 1e-04     0 1.0000000
#> 46:   Inf  0.98193889     0  0.99237531  0.99237531 1e-04     0 0.9470593
#> 47:   Inf  0.99475000     0  0.99570370  0.99570370 1e-04     0 1.0000000
#> 48:   Inf  0.98178889     0  0.99102160  0.99102160 1e-04     0 0.5022123
#>     GEC50      GR Max  ec50      GR Inf        GR 0  h GR E Inf   p_value
#>     <num>       <num> <num>       <num>       <num> <num> <num>     <num>
#>              rss Concentration N_conc         E0             fit_type
#>            <num>         <num>  <int>      <num>               <char>
#>  1: 0.0008866115             1      9 0.01649918 DRCConstantFitResult
#>  2: 0.0004315725             1      9 0.01839689 DRCConstantFitResult
#>  3: 0.0007135742             1      9 0.02138625 DRCConstantFitResult
#>  4: 0.0012209428             1      9 0.02595851 DRCConstantFitResult
#>  5: 0.0022302917             1      9 0.02747480 DRCConstantFitResult
#>  6: 0.0113145908             1      9 0.03277132 DRCConstantFitResult
#>  7: 0.0033764235             1      9 0.02448461 DRCConstantFitResult
#>  8: 0.0078975847             1      9 0.03212551 DRCConstantFitResult
#>  9: 0.0007231072             1      9 0.02728644 DRCConstantFitResult
#> 10: 0.0027121533             1      9 0.02315973 DRCConstantFitResult
#> 11: 0.0016184375             1      9 0.01779364 DRCConstantFitResult
#> 12: 0.0018844868             1      9 0.02788363 DRCConstantFitResult
#> 13: 0.0034670347             1      9 0.04144568 DRCConstantFitResult
#> 14: 0.0063340983             1      9 0.01543815 DRCConstantFitResult
#> 15: 0.0014072769             1      9 0.01812046 DRCConstantFitResult
#> 16: 0.0007678705             1      9 0.02714955 DRCConstantFitResult
#> 17: 0.0012846424             1      9 0.01975923 DRCConstantFitResult
#> 18: 0.0027127052             1      9 0.03229288 DRCConstantFitResult
#> 19: 0.0009149403             1      9 0.01314100 DRCConstantFitResult
#> 20: 0.0003673962             1      9 0.01333922 DRCConstantFitResult
#> 21: 0.0001904803             1      9 0.01667500 DRCConstantFitResult
#> 22: 0.0017873127             1      9 0.02152308 DRCConstantFitResult
#> 23: 0.0015366333             1      9 0.02457477 DRCConstantFitResult
#> 24: 0.0009474410             1      9 0.01680536 DRCConstantFitResult
#> 25: 0.0001913761             1      9 0.02746480 DRCConstantFitResult
#> 26: 0.0003942057             1      9 0.02873235 DRCConstantFitResult
#> 27: 0.0002641617             1      9 0.02824695 DRCConstantFitResult
#> 28: 0.0001847002             1      9 0.02798887 DRCConstantFitResult
#> 29: 0.0003173711             1      9 0.02860886 DRCConstantFitResult
#> 30: 0.0004066196             1      9 0.02893784 DRCConstantFitResult
#> 31: 0.0002318267             1      9 0.02974214 DRCConstantFitResult
#> 32: 0.0006522847             1      9 0.02977052 DRCConstantFitResult
#> 33: 0.0002597482             1      9 0.02978818 DRCConstantFitResult
#> 34: 0.0009760935             1      9 0.02804300 DRCConstantFitResult
#> 35: 0.0003425882             1      9 0.02895447 DRCConstantFitResult
#> 36: 0.0001608500             1      9 0.02855131 DRCConstantFitResult
#> 37: 0.0003994238             1      9 0.02823440 DRCConstantFitResult
#> 38: 0.0003875042             1      9 0.02856573 DRCConstantFitResult
#> 39: 0.0009228802             1      9 0.02889972 DRCConstantFitResult
#> 40: 0.0004802806             1      9 0.03010221 DRCConstantFitResult
#> 41: 0.0002059010             1      9 0.03105234 DRCConstantFitResult
#> 42: 0.0002710891             1      9 0.02977945 DRCConstantFitResult
#> 43: 0.0007010440             1      9 0.02881507 DRCConstantFitResult
#> 44: 0.0004884894             1      9 0.02909892 DRCConstantFitResult
#> 45: 0.0005460898             1      9 0.02845883 DRCConstantFitResult
#> 46: 0.0005124332             1      9 0.02799662 DRCConstantFitResult
#> 47: 0.0004795990             1      9 0.02781475 DRCConstantFitResult
#> 48: 0.0002233779             1      9 0.02960396 DRCConstantFitResult
#>              rss Concentration N_conc         E0             fit_type
#>            <num>         <num>  <int>      <num>               <char>
#>     normalization_type fit_source cotrt_value ratio dilution_drug Gnumber
#>                 <char>     <char>       <num> <num>        <char>  <char>
#>  1:                 RV        gDR       3.160    NA        drug_2  G00001
#>  2:                 RV        gDR       3.160    NA        drug_2  G00002
#>  3:                 RV        gDR      10.000    NA        drug_2  G00011
#>  4:                 GR        gDR       1.000    NA        drug_2  G00001
#>  5:                 GR        gDR       0.316    NA        drug_2  G00002
#>  6:                 GR        gDR       3.160    NA        drug_2  G00011
#>  7:                 GR        gDR      10.000    NA        drug_2  G00001
#>  8:                 GR        gDR       3.160    NA        drug_2  G00002
#>  9:                 GR        gDR       3.160    NA        drug_2  G00011
#> 10:                 GR        gDR      10.000    NA        drug_2  G00001
#> 11:                 GR        gDR       0.316    NA        drug_2  G00002
#> 12:                 GR        gDR       1.000    NA        drug_2  G00011
#> 13:                 GR        gDR      10.000    NA        drug_2  G00001
#> 14:                 GR        gDR      10.000    NA        drug_2  G00002
#> 15:                 GR        gDR      10.000    NA        drug_2  G00011
#> 16:                 GR        gDR      10.000    NA        drug_2  G00001
#> 17:                 GR        gDR      10.000    NA        drug_2  G00002
#> 18:                 GR        gDR      10.000    NA        drug_2  G00011
#> 19:                 GR        gDR      10.000    NA        drug_2  G00001
#> 20:                 GR        gDR      10.000    NA        drug_2  G00002
#> 21:                 GR        gDR      10.000    NA        drug_2  G00011
#> 22:                 GR        gDR       3.160    NA        drug_2  G00001
#> 23:                 GR        gDR       3.160    NA        drug_2  G00002
#> 24:                 GR        gDR      10.000    NA        drug_2  G00011
#> 25:                 RV        gDR          NA    NA          <NA>  G00021
#> 26:                 RV        gDR          NA    NA          <NA>  G00026
#> 27:                 RV        gDR          NA    NA          <NA>  G00031
#> 28:                 RV        gDR          NA    NA          <NA>  G00021
#> 29:                 RV        gDR          NA    NA          <NA>  G00026
#> 30:                 RV        gDR          NA    NA          <NA>  G00031
#> 31:                 RV        gDR          NA    NA          <NA>  G00021
#> 32:                 RV        gDR          NA    NA          <NA>  G00026
#> 33:                 RV        gDR          NA    NA          <NA>  G00031
#> 34:                 RV        gDR          NA    NA          <NA>  G00021
#> 35:                 RV        gDR          NA    NA          <NA>  G00026
#> 36:                 RV        gDR          NA    NA          <NA>  G00031
#> 37:                 RV        gDR          NA    NA          <NA>  G00021
#> 38:                 RV        gDR          NA    NA          <NA>  G00026
#> 39:                 RV        gDR          NA    NA          <NA>  G00031
#> 40:                 RV        gDR          NA    NA          <NA>  G00021
#> 41:                 RV        gDR          NA    NA          <NA>  G00026
#> 42:                 RV        gDR          NA    NA          <NA>  G00031
#> 43:                 RV        gDR          NA    NA          <NA>  G00021
#> 44:                 RV        gDR          NA    NA          <NA>  G00026
#> 45:                 RV        gDR          NA    NA          <NA>  G00031
#> 46:                 RV        gDR          NA    NA          <NA>  G00021
#> 47:                 RV        gDR          NA    NA          <NA>  G00026
#> 48:                 RV        gDR          NA    NA          <NA>  G00031
#>     normalization_type fit_source cotrt_value ratio dilution_drug Gnumber
#>                 <char>     <char>       <num> <num>        <char>  <char>
#>     Drug MOA Gnumber_2 Drug Name 2 Drug MOA2 Duration    clid   Tissue
#>       <char>    <char>      <char>    <char>    <num>  <char>   <char>
#>  1:    moa_A    G00021    drug_021     moa_D       72 CL00010 tissue_x
#>  2:    moa_A    G00021    drug_021     moa_D       72 CL00010 tissue_x
#>  3:    moa_B    G00021    drug_021     moa_D       72 CL00010 tissue_x
#>  4:    moa_A    G00021    drug_021     moa_D       72 CL00014 tissue_x
#>  5:    moa_A    G00021    drug_021     moa_D       72 CL00014 tissue_x
#>  6:    moa_B    G00021    drug_021     moa_D       72 CL00014 tissue_x
#>  7:    moa_A    G00021    drug_021     moa_D       72 CL00018 tissue_y
#>  8:    moa_A    G00021    drug_021     moa_D       72 CL00018 tissue_y
#>  9:    moa_B    G00021    drug_021     moa_D       72 CL00018 tissue_y
#> 10:    moa_A    G00021    drug_021     moa_D       72 CL00022 tissue_z
#> 11:    moa_A    G00021    drug_021     moa_D       72 CL00022 tissue_z
#> 12:    moa_B    G00021    drug_021     moa_D       72 CL00022 tissue_z
#> 13:    moa_A    G00021    drug_021     moa_D       72 CL00026 tissue_w
#> 14:    moa_A    G00021    drug_021     moa_D       72 CL00026 tissue_w
#> 15:    moa_B    G00021    drug_021     moa_D       72 CL00026 tissue_w
#> 16:    moa_A    G00021    drug_021     moa_D       72 CL00030 tissue_w
#> 17:    moa_A    G00021    drug_021     moa_D       72 CL00030 tissue_w
#> 18:    moa_B    G00021    drug_021     moa_D       72 CL00030 tissue_w
#> 19:    moa_A    G00021    drug_021     moa_D       72 CL00034 tissue_w
#> 20:    moa_A    G00021    drug_021     moa_D       72 CL00034 tissue_w
#> 21:    moa_B    G00021    drug_021     moa_D       72 CL00034 tissue_w
#> 22:    moa_A    G00021    drug_021     moa_D       72 CL00038 tissue_w
#> 23:    moa_A    G00021    drug_021     moa_D       72 CL00038 tissue_w
#> 24:    moa_B    G00021    drug_021     moa_D       72 CL00038 tissue_w
#> 25:    moa_D      <NA>        <NA>      <NA>       72 CL00010 tissue_x
#> 26:    moa_E      <NA>        <NA>      <NA>       72 CL00010 tissue_x
#> 27:    moa_F      <NA>        <NA>      <NA>       72 CL00010 tissue_x
#> 28:    moa_D      <NA>        <NA>      <NA>       72 CL00014 tissue_x
#> 29:    moa_E      <NA>        <NA>      <NA>       72 CL00014 tissue_x
#> 30:    moa_F      <NA>        <NA>      <NA>       72 CL00014 tissue_x
#> 31:    moa_D      <NA>        <NA>      <NA>       72 CL00018 tissue_y
#> 32:    moa_E      <NA>        <NA>      <NA>       72 CL00018 tissue_y
#> 33:    moa_F      <NA>        <NA>      <NA>       72 CL00018 tissue_y
#> 34:    moa_D      <NA>        <NA>      <NA>       72 CL00022 tissue_z
#> 35:    moa_E      <NA>        <NA>      <NA>       72 CL00022 tissue_z
#> 36:    moa_F      <NA>        <NA>      <NA>       72 CL00022 tissue_z
#> 37:    moa_D      <NA>        <NA>      <NA>       72 CL00026 tissue_w
#> 38:    moa_E      <NA>        <NA>      <NA>       72 CL00026 tissue_w
#> 39:    moa_F      <NA>        <NA>      <NA>       72 CL00026 tissue_w
#> 40:    moa_D      <NA>        <NA>      <NA>       72 CL00030 tissue_w
#> 41:    moa_E      <NA>        <NA>      <NA>       72 CL00030 tissue_w
#> 42:    moa_F      <NA>        <NA>      <NA>       72 CL00030 tissue_w
#> 43:    moa_D      <NA>        <NA>      <NA>       72 CL00034 tissue_w
#> 44:    moa_E      <NA>        <NA>      <NA>       72 CL00034 tissue_w
#> 45:    moa_F      <NA>        <NA>      <NA>       72 CL00034 tissue_w
#> 46:    moa_D      <NA>        <NA>      <NA>       72 CL00038 tissue_w
#> 47:    moa_E      <NA>        <NA>      <NA>       72 CL00038 tissue_w
#> 48:    moa_F      <NA>        <NA>      <NA>       72 CL00038 tissue_w
#>     Drug MOA Gnumber_2 Drug Name 2 Drug MOA2 Duration    clid   Tissue
#>       <char>    <char>      <char>    <char>    <num>  <char>   <char>
#>     ReferenceDivisionTime        GR50        IC50       E Max    GR value
#>                     <num>       <num>       <num>       <num>       <num>
#>  1:                    22  0.31359259  0.31359259  0.31359259  0.31359259
#>  2:                    22  0.46488148  0.46488148  0.46488148  0.46488148
#>  3:                    22  0.30005185  0.30005185  0.30005185  0.30005185
#>  4:                    38  0.54018519  0.54018519  0.54018519  0.54018519
#>  5:                    38  0.46981481  0.46981481  0.46981481  0.46981481
#>  6:                    38 -0.27716296 -0.27716296 -0.27716296 -0.27716296
#>  7:                    54 -0.70610370 -0.70610370 -0.70610370 -0.70610370
#>  8:                    54 -0.86958519 -0.86958519 -0.86958519 -0.86958519
#>  9:                    54 -0.60236296 -0.60236296 -0.60236296 -0.60236296
#> 10:                    70 -0.73374444 -0.73374444 -0.73374444 -0.73374444
#> 11:                    70 -0.58530000 -0.58530000 -0.58530000 -0.58530000
#> 12:                    70 -0.54707037 -0.54707037 -0.54707037 -0.54707037
#> 13:                    26 -0.27410000 -0.27410000 -0.27410000 -0.27410000
#> 14:                    26  0.19546296  0.19546296  0.19546296  0.19546296
#> 15:                    26  0.08634815  0.08634815  0.08634815  0.08634815
#> 16:                    42 -0.36797407 -0.36797407 -0.36797407 -0.36797407
#> 17:                    42 -0.30042222 -0.30042222 -0.30042222 -0.30042222
#> 18:                    42 -0.41170741 -0.41170741 -0.41170741 -0.41170741
#> 19:                    58 -0.69018148 -0.69018148 -0.69018148 -0.69018148
#> 20:                    58 -0.64542222 -0.64542222 -0.64542222 -0.64542222
#> 21:                    58 -0.64786296 -0.64786296 -0.64786296 -0.64786296
#> 22:                    74 -0.75242963 -0.75242963 -0.75242963 -0.75242963
#> 23:                    74 -0.54112963 -0.54112963 -0.54112963 -0.54112963
#> 24:                    74 -0.43574815 -0.43574815 -0.43574815 -0.43574815
#> 25:                    22  0.99633210  0.99633210  0.99633210  0.99633210
#> 26:                    22  0.99605432  0.99605432  0.99605432  0.99605432
#> 27:                    22  0.99742531  0.99742531  0.99742531  0.99742531
#> 28:                    38  0.99882160  0.99882160  0.99882160  0.99882160
#> 29:                    38  1.00082469  1.00082469  1.00082469  1.00082469
#> 30:                    38  0.99524444  0.99524444  0.99524444  0.99524444
#> 31:                    54  0.99930432  0.99930432  0.99930432  0.99930432
#> 32:                    54  0.99741975  0.99741975  0.99741975  0.99741975
#> 33:                    54  0.99455556  0.99455556  0.99455556  0.99455556
#> 34:                    70  0.99650309  0.99650309  0.99650309  0.99650309
#> 35:                    70  1.00027901  1.00027901  1.00027901  1.00027901
#> 36:                    70  0.99801481  0.99801481  0.99801481  0.99801481
#> 37:                    26  0.99741235  0.99741235  0.99741235  0.99741235
#> 38:                    26  1.00112284  1.00112284  1.00112284  1.00112284
#> 39:                    26  0.99772222  0.99772222  0.99772222  0.99772222
#> 40:                    42  0.99697037  0.99697037  0.99697037  0.99697037
#> 41:                    42  1.00156420  1.00156420  1.00156420  1.00156420
#> 42:                    42  0.99700185  0.99700185  0.99700185  0.99700185
#> 43:                    58  1.00209074  1.00209074  1.00209074  1.00209074
#> 44:                    58  1.00434630  1.00434630  1.00434630  1.00434630
#> 45:                    58  1.00456852  1.00456852  1.00456852  1.00456852
#> 46:                    74  0.99237531  0.99237531  0.99237531  0.99237531
#> 47:                    74  0.99570370  0.99570370  0.99570370  0.99570370
#> 48:                    74  0.99102160  0.99102160  0.99102160  0.99102160
#>     ReferenceDivisionTime        GR50        IC50       E Max    GR value
#>                     <num>       <num>       <num>       <num>       <num>
#> 
#> $raw_data
#>                                                  rId
#>                                               <char>
#>    1: G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#>    2: G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#>    3: G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#>    4: G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#>    5: G00001_drug_001_moa_A_G00021_drug_021_moa_D_72
#>   ---                                               
#> 4268:                       G00021_drug_021_moa_D_72
#> 4269:                       G00026_drug_026_moa_E_72
#> 4270:                       G00026_drug_026_moa_E_72
#> 4271:                       G00031_drug_031_moa_F_72
#> 4272:                       G00031_drug_031_moa_F_72
#>                                   cId    x_mean       x_AOC x_AOC_range  xc50
#>                                <char>     <num>       <num>       <num> <num>
#>    1: CL00010_cellline_AA_tissue_x_22 0.3135926 0.686407407 0.686407407  -Inf
#>    2: CL00010_cellline_AA_tissue_x_22 0.3142704 0.685729630 0.685729630  -Inf
#>    3: CL00010_cellline_AA_tissue_x_22 0.3148000 0.685200000 0.685200000  -Inf
#>    4: CL00010_cellline_AA_tissue_x_22 0.3158481 0.684151852 0.684151852  -Inf
#>    5: CL00010_cellline_AA_tissue_x_22 0.3323388 0.667661219 0.601409743  -Inf
#>   ---                                                                        
#> 4268: CL00038_cellline_NE_tissue_w_74 0.9843549 0.015645062 0.015645062   Inf
#> 4269: CL00038_cellline_NE_tissue_w_74 0.9957037 0.004296296 0.004296296   Inf
#> 4270: CL00038_cellline_NE_tissue_w_74 0.9911914 0.008808642 0.008808642   Inf
#> 4271: CL00038_cellline_NE_tissue_w_74 0.9910216 0.008978395 0.008978395   Inf
#> 4272: CL00038_cellline_NE_tissue_w_74 0.9815704 0.018429630 0.018429630   Inf
#>           x_max   ec50     x_inf       x_0         h        r2    p_value
#>           <num>  <num>     <num>     <num>     <num>     <num>      <num>
#>    1: 0.3110667 0.0000 0.3135926 0.3135926 0.0001000  0.000000 1.00000000
#>    2: 0.2997000 0.0000 0.3142704 0.3142704 0.0001000  0.000000 0.09012268
#>    3: 0.3055667 0.0000 0.3148000 0.3148000 0.0001000  0.000000 0.48014356
#>    4: 0.3068667 0.0000 0.3158481 0.3158481 0.0001000  0.000000 0.59634451
#>    5: 0.2999333 0.0101 0.0000000 1.0000000 0.1530207 -9.123791 1.00000000
#>   ---                                                                    
#> 4268: 0.9629056 0.0000 0.9843549 0.9843549 0.0001000  0.000000 0.95421110
#> 4269: 0.9947500 0.0000 0.9957037 0.9957037 0.0001000  0.000000 1.00000000
#> 4270: 0.9892333 0.0000 0.9911914 0.9911914 0.0001000  0.000000 1.00000000
#> 4271: 0.9817889 0.0000 0.9910216 0.9910216 0.0001000  0.000000 0.50221225
#> 4272: 0.9625889 0.0000 0.9815704 0.9815704 0.0001000  0.000000 0.49077639
#>                rss maxlog10Concentration N_conc   x_sd_avg
#>              <num>                 <num>  <int>      <num>
#>    1: 0.0008866115              1.000000      9 0.01649918
#>    2: 0.0003075320              1.000000      9 0.02110240
#>    3: 0.0003981517              1.000000      9 0.01675435
#>    4: 0.0008042315              1.000000      9 0.01691745
#>    5: 0.0117634036              1.004321      5 0.01551077
#>   ---                                                     
#> 4268: 0.0021689801              1.000000      9 0.05753726
#> 4269: 0.0004795990              1.000000      9 0.02781475
#> 4270: 0.0021142487              1.000000      9 0.05716746
#> 4271: 0.0002233779              1.000000      9 0.02960396
#> 4272: 0.0009375695              1.000000      9 0.06084279
#>                     fit_type normalization_type fit_source cotrt_value ratio
#>                       <char>             <char>     <char>       <num> <num>
#>    1:   DRCConstantFitResult                 RV        gDR       3.160    NA
#>    2:   DRCConstantFitResult                 RV        gDR      10.000    NA
#>    3:   DRCConstantFitResult                 RV        gDR       0.316    NA
#>    4:   DRCConstantFitResult                 RV        gDR       1.000    NA
#>    5: DRC3pHillFitModelFixS0                 RV        gDR          NA  0.01
#>   ---                                                                       
#> 4268:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 4269:   DRCConstantFitResult                 RV        gDR          NA    NA
#> 4270:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 4271:   DRCConstantFitResult                 RV        gDR          NA    NA
#> 4272:   DRCConstantFitResult                 GR        gDR          NA    NA
#>       dilution_drug Gnumber DrugName drug_moa Gnumber_2 DrugName_2 drug_moa_2
#>              <char>  <char>   <char>   <char>    <char>     <char>     <char>
#>    1:        drug_2  G00001 drug_001    moa_A    G00021   drug_021      moa_D
#>    2:        drug_2  G00001 drug_001    moa_A    G00021   drug_021      moa_D
#>    3:        drug_2  G00001 drug_001    moa_A    G00021   drug_021      moa_D
#>    4:        drug_2  G00001 drug_001    moa_A    G00021   drug_021      moa_D
#>    5:    codilution  G00001 drug_001    moa_A    G00021   drug_021      moa_D
#>   ---                                                                        
#> 4268:          <NA>  G00021 drug_021    moa_D      <NA>       <NA>       <NA>
#> 4269:          <NA>  G00026 drug_026    moa_E      <NA>       <NA>       <NA>
#> 4270:          <NA>  G00026 drug_026    moa_E      <NA>       <NA>       <NA>
#> 4271:          <NA>  G00031 drug_031    moa_F      <NA>       <NA>       <NA>
#> 4272:          <NA>  G00031 drug_031    moa_F      <NA>       <NA>       <NA>
#>       Duration    clid CellLineName   Tissue ReferenceDivisionTime
#>          <num>  <char>       <char>   <char>                 <num>
#>    1:       72 CL00010  cellline_AA tissue_x                    22
#>    2:       72 CL00010  cellline_AA tissue_x                    22
#>    3:       72 CL00010  cellline_AA tissue_x                    22
#>    4:       72 CL00010  cellline_AA tissue_x                    22
#>    5:       72 CL00010  cellline_AA tissue_x                    22
#>   ---                                                             
#> 4268:       72 CL00038  cellline_NE tissue_w                    74
#> 4269:       72 CL00038  cellline_NE tissue_w                    74
#> 4270:       72 CL00038  cellline_NE tissue_w                    74
#> 4271:       72 CL00038  cellline_NE tissue_w                    74
#> 4272:       72 CL00038  cellline_NE tissue_w                    74
#> 
```
