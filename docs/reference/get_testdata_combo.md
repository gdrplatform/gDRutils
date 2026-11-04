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
#>  1: CL00010_cellline_AA_tissue_x_22  0.53767650  0.4623234986  0.5147181921
#>  2: CL00010_cellline_AA_tissue_x_22  0.72322487  0.2767751274  0.3197659339
#>  3: CL00010_cellline_AA_tissue_x_22  0.53746769  0.4625323137  0.5320801477
#>  4: CL00014_cellline_EA_tissue_x_38  0.69147901  0.3085209866  0.3631789686
#>  5: CL00014_cellline_EA_tissue_x_38  0.64379021  0.3562097912  0.4199806960
#>  6: CL00014_cellline_EA_tissue_x_38  0.27739831  0.7226016893  0.8280503313
#>  7: CL00018_cellline_IB_tissue_y_54 -0.08545482  1.0854548199  1.2681759314
#>  8: CL00018_cellline_IB_tissue_y_54 -0.54773187  1.5477318734  1.7511370571
#>  9: CL00018_cellline_IB_tissue_y_54  0.09674213  0.9032578664  1.0327415031
#> 10: CL00022_cellline_MC_tissue_z_70 -0.04151376  1.0415137626  1.2105084802
#> 11: CL00022_cellline_MC_tissue_z_70 -0.09405214  1.0940521444  1.2926430848
#> 12: CL00022_cellline_MC_tissue_z_70 -0.02599303  1.0259930262  1.2086461204
#> 13: CL00026_cellline_BC_tissue_w_26  0.37086767  0.6291323319  0.7071127226
#> 14: CL00026_cellline_BC_tissue_w_26  0.57264538  0.4273546175  0.4831064806
#> 15: CL00026_cellline_BC_tissue_w_26  0.58639343  0.4136065672  0.4603668126
#> 16: CL00030_cellline_FD_tissue_w_42  0.09685607  0.9031439304  1.0501402230
#> 17: CL00030_cellline_FD_tissue_w_42  0.14894011  0.8510598861  0.9898112854
#> 18: CL00030_cellline_FD_tissue_w_42  0.24635719  0.7536428124  0.8626105539
#> 19: CL00034_cellline_JE_tissue_w_58  0.28694710  0.7130528968  0.7805937021
#> 20: CL00034_cellline_JE_tissue_w_58  0.34689033  0.6531096718  0.7063641569
#> 21: CL00034_cellline_JE_tissue_w_58  0.37678374  0.6232162627  0.6648129271
#> 22: CL00038_cellline_NE_tissue_w_74 -0.11338165  1.1133816528  1.3003526511
#> 23: CL00038_cellline_NE_tissue_w_74  0.10270596  0.8972940359  1.0366366127
#> 24: CL00038_cellline_NE_tissue_w_74  0.27217808  0.7278219152  0.8220826277
#> 25: CL00010_cellline_AA_tissue_x_22  0.99759568  0.0024043210  0.0024043210
#> 26: CL00010_cellline_AA_tissue_x_22  0.99741111  0.0025888889  0.0025888889
#> 27: CL00010_cellline_AA_tissue_x_22  0.99825741  0.0017425926  0.0017425926
#> 28: CL00014_cellline_EA_tissue_x_38  0.99856358  0.0014364198  0.0014364198
#> 29: CL00014_cellline_EA_tissue_x_38  1.00066358 -0.0006635802 -0.0006635802
#> 30: CL00014_cellline_EA_tissue_x_38  0.99476173  0.0052382716  0.0052382716
#> 31: CL00018_cellline_IB_tissue_y_54  0.99878642  0.0012135802  0.0012135802
#> 32: CL00018_cellline_IB_tissue_y_54  0.99595123  0.0040487654  0.0040487654
#> 33: CL00018_cellline_IB_tissue_y_54  0.99165370  0.0083462963  0.0083462963
#> 34: CL00022_cellline_MC_tissue_z_70  0.99317531  0.0068246914  0.0068246914
#> 35: CL00022_cellline_MC_tissue_z_70  1.00052407 -0.0005240741 -0.0005240741
#> 36: CL00022_cellline_MC_tissue_z_70  0.99611605  0.0038839506  0.0038839506
#> 37: CL00026_cellline_BC_tissue_w_26  0.99794321  0.0020567901  0.0020567901
#> 38: CL00026_cellline_BC_tissue_w_26  1.00062222 -0.0006222222 -0.0006222222
#> 39: CL00026_cellline_BC_tissue_w_26  0.99814815  0.0018518519  0.0018518519
#> 40: CL00030_cellline_FD_tissue_w_42  0.99623827  0.0037617284  0.0037617284
#> 41: CL00030_cellline_FD_tissue_w_42  1.00159568 -0.0015956790 -0.0015956790
#> 42: CL00030_cellline_FD_tissue_w_42  0.99628765  0.0037123457  0.0037123457
#> 43: CL00034_cellline_JE_tissue_w_58  1.00323025 -0.0032302469 -0.0032302469
#> 44: CL00034_cellline_JE_tissue_w_58  1.00686605 -0.0068660494 -0.0068660494
#> 45: CL00034_cellline_JE_tissue_w_58  1.00722840 -0.0072283951 -0.0072283951
#> 46: CL00038_cellline_NE_tissue_w_74  0.98435494  0.0156450617  0.0156450617
#> 47: CL00038_cellline_NE_tissue_w_74  0.99119136  0.0088086420  0.0088086420
#> 48: CL00038_cellline_NE_tissue_w_74  0.98157037  0.0184296296  0.0184296296
#>                                 cId        EC50         x_AOC          h RV
#>                              <char>       <num>         <num>         <num>
#>           GEC50      GR Max        ec50     GR Inf      GR 0      h GR
#>           <num>       <num>       <num>      <num>     <num>     <num>
#>  1: 0.041757091  0.39526667 0.005734645  0.4106514 0.9975957 0.8649596
#>  2:         Inf  0.57806667 0.024255371  0.5757585 0.9975957 1.1872553
#>  3: 0.030946136  0.35043333 0.009614897  0.3814293 0.9975957 1.2269965
#>  4:         Inf  0.52480000 0.022049220  0.5350320 0.9985636 2.2737618
#>  5: 0.095719797  0.48110000 0.016855329  0.4853263 0.9985636 2.0300315
#>  6: 0.044862453 -0.34343333 0.062594800 -0.3120684 0.9985636 1.4646692
#>  7: 0.012042965 -0.69040000 0.020661474 -0.6223416 0.9987864 1.5024257
#>  8: 0.002133997 -0.89966667 0.003772607 -0.8277781 0.9987864 1.7183863
#>  9: 0.026057418 -0.60593333 0.052573808 -0.5885658 0.9987864 1.1118692
#> 10: 0.023123440 -0.71603333 0.037959719 -0.7197169 0.9931753 1.8268068
#> 11: 0.012230000 -0.58916667 0.017075863 -0.5826234 0.9931753 2.3557192
#> 12: 0.016238640 -0.54576667 0.021870787 -0.5431426 0.9931753 2.5158950
#> 13: 0.090456250 -0.31593333 0.116104772 -0.2977933 0.9979432 1.8882736
#> 14: 0.125063237  0.17123333 0.085432582  0.1749013 0.9979432 1.1187779
#> 15: 0.160534855  0.09876667 0.145431058  0.1027648 0.9979432 2.2868119
#> 16: 0.016162342 -0.39170000 0.024713623 -0.3887336 0.9962383 1.3722194
#> 17: 0.017780566 -0.30823333 0.025377048 -0.3140850 0.9962383 1.3914933
#> 18: 0.055729164 -0.40870000 0.070539244 -0.3988467 0.9962383 2.5207591
#> 19: 0.128498043 -0.69866667 0.206313751 -0.6933637 1.0032302 1.8236788
#> 20: 0.169904529 -0.63663333 0.253292959 -0.6371314 1.0032302 2.0415618
#> 21: 0.200757187 -0.64683333 0.313087466 -0.6567172 1.0032302 1.8729047
#> 22: 0.015800759 -0.74076667 0.026711978 -0.7262111 0.9843549 1.7690867
#> 23: 0.031860530 -0.54476667 0.047066519 -0.5331529 0.9843549 1.9414767
#> 24: 0.065854826 -0.43410000 0.095929807 -0.4273946 0.9843549 1.7268028
#> 25:         Inf  0.99263889 0.000000000  0.9975957 1.0000000 0.0001000
#> 26:         Inf  0.99342222 0.000000000  0.9974111 1.0000000 0.0001000
#> 27:         Inf  0.99976667 0.000000000  0.9982574 1.0000000 0.0001000
#> 28:         Inf  1.00167222 0.000000000  0.9985636 1.0000000 0.0001000
#> 29:         Inf  0.98942222 0.000000000  1.0006636 1.0000000 0.0001000
#> 30:         Inf  0.98533333 0.000000000  0.9947617 1.0000000 0.0001000
#> 31:         Inf  0.99144444 0.000000000  0.9987864 1.0000000 0.0001000
#> 32:         Inf  0.97545556 0.000000000  0.9959512 1.0000000 0.0001000
#> 33:         Inf  0.98621667 0.000000000  0.9916537 1.0000000 0.0001000
#> 34:         Inf  0.98348889 0.000000000  0.9931753 1.0000000 0.0001000
#> 35:         Inf  1.00531667 0.000000000  1.0005241 1.0000000 0.0001000
#> 36:         Inf  0.99403889 0.000000000  0.9961160 1.0000000 0.0001000
#> 37:         Inf  0.99560000 0.000000000  0.9979432 1.0000000 0.0001000
#> 38:         Inf  0.99286667 0.000000000  1.0006222 1.0000000 0.0001000
#> 39:         Inf  0.98954444 0.000000000  0.9981481 1.0000000 0.0001000
#> 40:         Inf  0.99541111 0.000000000  0.9962383 1.0000000 0.0001000
#> 41:         Inf  0.99352222 0.000000000  1.0015957 1.0000000 0.0001000
#> 42:         Inf  0.99225000 0.000000000  0.9962877 1.0000000 0.0001000
#> 43:         Inf  0.99173333 0.000000000  1.0032302 1.0000000 0.0001000
#> 44:         Inf  0.99821111 0.000000000  1.0068660 1.0000000 0.0001000
#> 45:         Inf  0.99041667 0.000000000  1.0072284 1.0000000 0.0001000
#> 46:         Inf  0.96290556 0.000000000  0.9843549 1.0000000 0.0001000
#> 47:         Inf  0.98923333 0.000000000  0.9911914 1.0000000 0.0001000
#> 48:         Inf  0.96258889 0.000000000  0.9815704 1.0000000 0.0001000
#>           GEC50      GR Max        ec50     GR Inf      GR 0      h GR
#>           <num>       <num>       <num>      <num>     <num>     <num>
#>         E Inf      p_value          rss Concentration N_conc         E0
#>         <num>        <num>        <num>         <num>  <int>      <num>
#>  1: 0.9726993 3.362070e-06 0.0075604110             1      9 0.02085838
#>  2: 0.9979975 3.593539e-10 0.0004975599             1      9 0.01489589
#>  3: 0.9825391 7.034472e-07 0.0080658121             1      9 0.03115615
#>  4: 0.9953919 6.642266e-09 0.0017011863             1      9 0.02234048
#>  5: 0.9967143 2.033298e-09 0.0013961805             1      9 0.02621187
#>  6: 0.9978480 4.623417e-10 0.0062411428             1      9 0.02996668
#>  7: 0.9883400 1.711765e-07 0.0473949600             1      9 0.04322901
#>  8: 0.9902549 9.136005e-08 0.0291332693             1      9 0.04932347
#>  9: 0.9906998 7.757386e-08 0.0353634981             1      9 0.03811409
#> 10: 0.9984051 1.620243e-10 0.0082091291             1      9 0.03329279
#> 11: 0.9983562 1.800990e-10 0.0068095430             1      9 0.02653348
#> 12: 0.9983227 1.932733e-10 0.0070511529             1      9 0.03067701
#> 13: 0.9974088 8.856059e-10 0.0078640779             1      9 0.01951018
#> 14: 0.9964917 2.557688e-09 0.0035070547             1      9 0.01512344
#> 15: 0.9995239 2.355069e-12 0.0007294436             1      9 0.01834673
#> 16: 0.9983373 1.874476e-10 0.0047931977             1      9 0.03302222
#> 17: 0.9965793 2.340925e-09 0.0088553531             1      9 0.03411108
#> 18: 0.9996246 1.024765e-12 0.0013828060             1      9 0.02608182
#> 19: 0.9994820 3.163126e-12 0.0026264052             1      9 0.02980156
#> 20: 0.9995910 1.383367e-12 0.0019526686             1      9 0.03882123
#> 21: 0.9998638 2.946370e-14 0.0006259533             1      9 0.03771691
#> 22: 0.9993618 6.567596e-12 0.0030586371             1      9 0.03563022
#> 23: 0.9984119 1.596252e-10 0.0067482025             1      9 0.02566952
#> 24: 0.9998139 8.791249e-14 0.0006463331             1      9 0.04093480
#> 25: 0.0000000           NA           NA             1      9 0.01683495
#> 26: 0.0000000           NA           NA             1      9 0.01759696
#> 27: 0.0000000           NA           NA             1      9 0.01728526
#> 28: 0.0000000           NA           NA             1      9 0.02956998
#> 29: 0.0000000           NA           NA             1      9 0.03023584
#> 30: 0.0000000           NA           NA             1      9 0.03059986
#> 31: 0.0000000           NA           NA             1      9 0.04462323
#> 32: 0.0000000           NA           NA             1      9 0.04469091
#> 33: 0.0000000           NA           NA             1      9 0.04471882
#> 34: 0.0000000           NA           NA             1      9 0.05453496
#> 35: 0.0000000           NA           NA             1      9 0.05629977
#> 36: 0.0000000           NA           NA             1      9 0.05551830
#> 37: 0.0000000           NA           NA             1      9 0.02040250
#> 38: 0.0000000           NA           NA             1      9 0.02063970
#> 39: 0.0000000           NA           NA             1      9 0.02092755
#> 40: 0.0000000           NA           NA             1      9 0.03513841
#> 41: 0.0000000           NA           NA             1      9 0.03621986
#> 42: 0.0000000           NA           NA             1      9 0.03476692
#> 43: 0.0000000           NA           NA             1      9 0.04638909
#> 44: 0.0000000           NA           NA             1      9 0.04684734
#> 45: 0.0000000           NA           NA             1      9 0.04582199
#> 46: 0.0000000           NA           NA             1      9 0.05753726
#> 47: 0.0000000           NA           NA             1      9 0.05716746
#> 48: 0.0000000           NA           NA             1      9 0.06084279
#>         E Inf      p_value          rss Concentration N_conc         E0
#>         <num>        <num>        <num>         <num>  <int>      <num>
#>                   fit_type normalization_type fit_source cotrt_value ratio
#>                     <char>             <char>     <char>       <num> <num>
#>  1: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#>  2: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#>  3: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#>  4: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#>  5: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#>  6: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#>  7: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#>  8: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#>  9: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 10: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 11: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 12: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 13: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 14: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 15: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 16: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 17: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 18: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 19: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 20: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 21: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 22: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 23: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 24: DRC3pHillFitModelFixS0                 GR        gDR       0.001    NA
#> 25:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 26:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 27:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 28:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 29:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 30:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 31:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 32:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 33:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 34:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 35:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 36:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 37:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 38:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 39:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 40:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 41:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 42:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 43:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 44:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 45:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 46:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 47:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 48:   DRCConstantFitResult                 GR        gDR          NA    NA
#>                   fit_type normalization_type fit_source cotrt_value ratio
#>                     <char>             <char>     <char>       <num> <num>
#>     dilution_drug Gnumber Drug MOA Gnumber_2 Drug Name 2 Drug MOA2 Duration
#>            <char>  <char>   <char>    <char>      <char>    <char>    <num>
#>  1:        drug_1  G00001    moa_A    G00021    drug_021     moa_D       72
#>  2:        drug_1  G00002    moa_A    G00021    drug_021     moa_D       72
#>  3:        drug_1  G00011    moa_B    G00021    drug_021     moa_D       72
#>  4:        drug_1  G00001    moa_A    G00021    drug_021     moa_D       72
#>  5:        drug_1  G00002    moa_A    G00021    drug_021     moa_D       72
#>  6:        drug_1  G00011    moa_B    G00021    drug_021     moa_D       72
#>  7:        drug_1  G00001    moa_A    G00021    drug_021     moa_D       72
#>  8:        drug_1  G00002    moa_A    G00021    drug_021     moa_D       72
#>  9:        drug_1  G00011    moa_B    G00021    drug_021     moa_D       72
#> 10:        drug_1  G00001    moa_A    G00021    drug_021     moa_D       72
#> 11:        drug_1  G00002    moa_A    G00021    drug_021     moa_D       72
#> 12:        drug_1  G00011    moa_B    G00021    drug_021     moa_D       72
#> 13:        drug_1  G00001    moa_A    G00021    drug_021     moa_D       72
#> 14:        drug_1  G00002    moa_A    G00021    drug_021     moa_D       72
#> 15:        drug_1  G00011    moa_B    G00021    drug_021     moa_D       72
#> 16:        drug_1  G00001    moa_A    G00021    drug_021     moa_D       72
#> 17:        drug_1  G00002    moa_A    G00021    drug_021     moa_D       72
#> 18:        drug_1  G00011    moa_B    G00021    drug_021     moa_D       72
#> 19:        drug_1  G00001    moa_A    G00021    drug_021     moa_D       72
#> 20:        drug_1  G00002    moa_A    G00021    drug_021     moa_D       72
#> 21:        drug_1  G00011    moa_B    G00021    drug_021     moa_D       72
#> 22:        drug_1  G00001    moa_A    G00021    drug_021     moa_D       72
#> 23:        drug_1  G00002    moa_A    G00021    drug_021     moa_D       72
#> 24:        drug_1  G00011    moa_B    G00021    drug_021     moa_D       72
#> 25:          <NA>  G00021    moa_D      <NA>        <NA>      <NA>       72
#> 26:          <NA>  G00026    moa_E      <NA>        <NA>      <NA>       72
#> 27:          <NA>  G00031    moa_F      <NA>        <NA>      <NA>       72
#> 28:          <NA>  G00021    moa_D      <NA>        <NA>      <NA>       72
#> 29:          <NA>  G00026    moa_E      <NA>        <NA>      <NA>       72
#> 30:          <NA>  G00031    moa_F      <NA>        <NA>      <NA>       72
#> 31:          <NA>  G00021    moa_D      <NA>        <NA>      <NA>       72
#> 32:          <NA>  G00026    moa_E      <NA>        <NA>      <NA>       72
#> 33:          <NA>  G00031    moa_F      <NA>        <NA>      <NA>       72
#> 34:          <NA>  G00021    moa_D      <NA>        <NA>      <NA>       72
#> 35:          <NA>  G00026    moa_E      <NA>        <NA>      <NA>       72
#> 36:          <NA>  G00031    moa_F      <NA>        <NA>      <NA>       72
#> 37:          <NA>  G00021    moa_D      <NA>        <NA>      <NA>       72
#> 38:          <NA>  G00026    moa_E      <NA>        <NA>      <NA>       72
#> 39:          <NA>  G00031    moa_F      <NA>        <NA>      <NA>       72
#> 40:          <NA>  G00021    moa_D      <NA>        <NA>      <NA>       72
#> 41:          <NA>  G00026    moa_E      <NA>        <NA>      <NA>       72
#> 42:          <NA>  G00031    moa_F      <NA>        <NA>      <NA>       72
#> 43:          <NA>  G00021    moa_D      <NA>        <NA>      <NA>       72
#> 44:          <NA>  G00026    moa_E      <NA>        <NA>      <NA>       72
#> 45:          <NA>  G00031    moa_F      <NA>        <NA>      <NA>       72
#> 46:          <NA>  G00021    moa_D      <NA>        <NA>      <NA>       72
#> 47:          <NA>  G00026    moa_E      <NA>        <NA>      <NA>       72
#> 48:          <NA>  G00031    moa_F      <NA>        <NA>      <NA>       72
#>     dilution_drug Gnumber Drug MOA Gnumber_2 Drug Name 2 Drug MOA2 Duration
#>            <char>  <char>   <char>    <char>      <char>    <char>    <num>
#>        clid   Tissue ReferenceDivisionTime        GR50        IC50       E Max
#>      <char>   <char>                 <num>       <num>       <num>       <num>
#>  1: CL00010 tissue_x                    22  0.53767650  0.53767650  0.53767650
#>  2: CL00010 tissue_x                    22  0.72322487  0.72322487  0.72322487
#>  3: CL00010 tissue_x                    22  0.53746769  0.53746769  0.53746769
#>  4: CL00014 tissue_x                    38  0.69147901  0.69147901  0.69147901
#>  5: CL00014 tissue_x                    38  0.64379021  0.64379021  0.64379021
#>  6: CL00014 tissue_x                    38  0.27739831  0.27739831  0.27739831
#>  7: CL00018 tissue_y                    54 -0.08545482 -0.08545482 -0.08545482
#>  8: CL00018 tissue_y                    54 -0.54773187 -0.54773187 -0.54773187
#>  9: CL00018 tissue_y                    54  0.09674213  0.09674213  0.09674213
#> 10: CL00022 tissue_z                    70 -0.04151376 -0.04151376 -0.04151376
#> 11: CL00022 tissue_z                    70 -0.09405214 -0.09405214 -0.09405214
#> 12: CL00022 tissue_z                    70 -0.02599303 -0.02599303 -0.02599303
#> 13: CL00026 tissue_w                    26  0.37086767  0.37086767  0.37086767
#> 14: CL00026 tissue_w                    26  0.57264538  0.57264538  0.57264538
#> 15: CL00026 tissue_w                    26  0.58639343  0.58639343  0.58639343
#> 16: CL00030 tissue_w                    42  0.09685607  0.09685607  0.09685607
#> 17: CL00030 tissue_w                    42  0.14894011  0.14894011  0.14894011
#> 18: CL00030 tissue_w                    42  0.24635719  0.24635719  0.24635719
#> 19: CL00034 tissue_w                    58  0.28694710  0.28694710  0.28694710
#> 20: CL00034 tissue_w                    58  0.34689033  0.34689033  0.34689033
#> 21: CL00034 tissue_w                    58  0.37678374  0.37678374  0.37678374
#> 22: CL00038 tissue_w                    74 -0.11338165 -0.11338165 -0.11338165
#> 23: CL00038 tissue_w                    74  0.10270596  0.10270596  0.10270596
#> 24: CL00038 tissue_w                    74  0.27217808  0.27217808  0.27217808
#> 25: CL00010 tissue_x                    22  0.99759568  0.99759568  0.99759568
#> 26: CL00010 tissue_x                    22  0.99741111  0.99741111  0.99741111
#> 27: CL00010 tissue_x                    22  0.99825741  0.99825741  0.99825741
#> 28: CL00014 tissue_x                    38  0.99856358  0.99856358  0.99856358
#> 29: CL00014 tissue_x                    38  1.00066358  1.00066358  1.00066358
#> 30: CL00014 tissue_x                    38  0.99476173  0.99476173  0.99476173
#> 31: CL00018 tissue_y                    54  0.99878642  0.99878642  0.99878642
#> 32: CL00018 tissue_y                    54  0.99595123  0.99595123  0.99595123
#> 33: CL00018 tissue_y                    54  0.99165370  0.99165370  0.99165370
#> 34: CL00022 tissue_z                    70  0.99317531  0.99317531  0.99317531
#> 35: CL00022 tissue_z                    70  1.00052407  1.00052407  1.00052407
#> 36: CL00022 tissue_z                    70  0.99611605  0.99611605  0.99611605
#> 37: CL00026 tissue_w                    26  0.99794321  0.99794321  0.99794321
#> 38: CL00026 tissue_w                    26  1.00062222  1.00062222  1.00062222
#> 39: CL00026 tissue_w                    26  0.99814815  0.99814815  0.99814815
#> 40: CL00030 tissue_w                    42  0.99623827  0.99623827  0.99623827
#> 41: CL00030 tissue_w                    42  1.00159568  1.00159568  1.00159568
#> 42: CL00030 tissue_w                    42  0.99628765  0.99628765  0.99628765
#> 43: CL00034 tissue_w                    58  1.00323025  1.00323025  1.00323025
#> 44: CL00034 tissue_w                    58  1.00686605  1.00686605  1.00686605
#> 45: CL00034 tissue_w                    58  1.00722840  1.00722840  1.00722840
#> 46: CL00038 tissue_w                    74  0.98435494  0.98435494  0.98435494
#> 47: CL00038 tissue_w                    74  0.99119136  0.99119136  0.99119136
#> 48: CL00038 tissue_w                    74  0.98157037  0.98157037  0.98157037
#>        clid   Tissue ReferenceDivisionTime        GR50        IC50       E Max
#>      <char>   <char>                 <num>       <num>       <num>       <num>
#>        GR value
#>           <num>
#>  1:  0.53767650
#>  2:  0.72322487
#>  3:  0.53746769
#>  4:  0.69147901
#>  5:  0.64379021
#>  6:  0.27739831
#>  7: -0.08545482
#>  8: -0.54773187
#>  9:  0.09674213
#> 10: -0.04151376
#> 11: -0.09405214
#> 12: -0.02599303
#> 13:  0.37086767
#> 14:  0.57264538
#> 15:  0.58639343
#> 16:  0.09685607
#> 17:  0.14894011
#> 18:  0.24635719
#> 19:  0.28694710
#> 20:  0.34689033
#> 21:  0.37678374
#> 22: -0.11338165
#> 23:  0.10270596
#> 24:  0.27217808
#> 25:  0.99759568
#> 26:  0.99741111
#> 27:  0.99825741
#> 28:  0.99856358
#> 29:  1.00066358
#> 30:  0.99476173
#> 31:  0.99878642
#> 32:  0.99595123
#> 33:  0.99165370
#> 34:  0.99317531
#> 35:  1.00052407
#> 36:  0.99611605
#> 37:  0.99794321
#> 38:  1.00062222
#> 39:  0.99814815
#> 40:  0.99623827
#> 41:  1.00159568
#> 42:  0.99628765
#> 43:  1.00323025
#> 44:  1.00686605
#> 45:  1.00722840
#> 46:  0.98435494
#> 47:  0.99119136
#> 48:  0.98157037
#>        GR value
#>           <num>
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
#>                                   cId    x_mean       x_AOC x_AOC_range
#>                                <char>     <num>       <num>       <num>
#>    1: CL00010_cellline_AA_tissue_x_22 0.5376765 0.462323499 0.514718192
#>    2: CL00010_cellline_AA_tissue_x_22 0.5235941 0.476405906 0.531279153
#>    3: CL00010_cellline_AA_tissue_x_22 0.5368776 0.463122395 0.513313373
#>    4: CL00010_cellline_AA_tissue_x_22 0.5297206 0.470279371 0.522997889
#>    5: CL00010_cellline_AA_tissue_x_22 0.5381584 0.461841592 0.515431114
#>   ---                                                                  
#> 4268: CL00038_cellline_NE_tissue_w_74 0.9923753 0.007624691 0.007624691
#> 4269: CL00038_cellline_NE_tissue_w_74 0.9911914 0.008808642 0.008808642
#> 4270: CL00038_cellline_NE_tissue_w_74 0.9957037 0.004296296 0.004296296
#> 4271: CL00038_cellline_NE_tissue_w_74 0.9815704 0.018429630 0.018429630
#> 4272: CL00038_cellline_NE_tissue_w_74 0.9910216 0.008978395 0.008978395
#>             xc50     x_max        ec50     x_inf       x_0         h        r2
#>            <num>     <num>       <num>     <num>     <num>     <num>     <num>
#>    1: 0.04175709 0.3952667 0.005734645 0.4106514 0.9975957 0.8649596 0.9726993
#>    2: 0.02908947 0.3972333 0.005100821 0.4042156 0.9975957 0.9464168 0.9840059
#>    3: 0.05011277 0.4060000 0.006456466 0.3958175 0.9975957 0.7630530 0.9747240
#>    4: 0.04087698 0.3842333 0.006325265 0.3910828 0.9975957 0.8141414 0.9679437
#>    5: 0.04217047 0.4137000 0.006222474 0.4060747 0.9975957 0.8713049 0.9767648
#>   ---                                                                         
#> 4268:        Inf 0.9819389 0.000000000 0.9923753 1.0000000 0.0001000 0.0000000
#> 4269:        Inf 0.9892333 0.000000000 0.9911914 1.0000000 0.0001000 0.0000000
#> 4270:        Inf 0.9947500 0.000000000 0.9957037 1.0000000 0.0001000 0.0000000
#> 4271:        Inf 0.9625889 0.000000000 0.9815704 1.0000000 0.0001000 0.0000000
#> 4272:        Inf 0.9817889 0.000000000 0.9910216 1.0000000 0.0001000 0.0000000
#>            p_value         rss maxlog10Concentration N_conc   x_sd_avg
#>              <num>       <num>                 <num>  <int>      <num>
#>    1: 3.362070e-06 0.007560411                     1      9 0.02085838
#>    2: 5.174404e-07 0.004483706                     1      9 0.01875524
#>    3: 2.567305e-06 0.006923634                     1      9 0.02656683
#>    4: 5.897883e-06 0.009512753                     1      9 0.01696155
#>    5: 1.912113e-06 0.006687606                     1      9 0.03165605
#>   ---                                                                 
#> 4268:           NA          NA                     1      9 0.02799662
#> 4269:           NA          NA                     1      9 0.05716746
#> 4270:           NA          NA                     1      9 0.02781475
#> 4271:           NA          NA                     1      9 0.06084279
#> 4272:           NA          NA                     1      9 0.02960396
#>                     fit_type normalization_type fit_source cotrt_value ratio
#>                       <char>             <char>     <char>       <num> <num>
#>    1: DRC3pHillFitModelFixS0                 GR        gDR     0.00100    NA
#>    2: DRC3pHillFitModelFixS0                 GR        gDR     0.00316    NA
#>    3: DRC3pHillFitModelFixS0                 GR        gDR     0.01000    NA
#>    4: DRC3pHillFitModelFixS0                 GR        gDR     0.03160    NA
#>    5: DRC3pHillFitModelFixS0                 GR        gDR     0.10000    NA
#>   ---                                                                       
#> 4268:   DRCConstantFitResult                 RV        gDR          NA    NA
#> 4269:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 4270:   DRCConstantFitResult                 RV        gDR          NA    NA
#> 4271:   DRCConstantFitResult                 GR        gDR          NA    NA
#> 4272:   DRCConstantFitResult                 RV        gDR          NA    NA
#>       dilution_drug Gnumber DrugName drug_moa Gnumber_2 DrugName_2 drug_moa_2
#>              <char>  <char>   <char>   <char>    <char>     <char>     <char>
#>    1:        drug_1  G00001 drug_001    moa_A    G00021   drug_021      moa_D
#>    2:        drug_1  G00001 drug_001    moa_A    G00021   drug_021      moa_D
#>    3:        drug_1  G00001 drug_001    moa_A    G00021   drug_021      moa_D
#>    4:        drug_1  G00001 drug_001    moa_A    G00021   drug_021      moa_D
#>    5:        drug_1  G00001 drug_001    moa_A    G00021   drug_021      moa_D
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
