# get_testdata_codilution

Function to obtain data from gDRtestData and prepare for unit tests

## Usage

``` r
get_testdata_codilution()
```

## Value

list with drugs, cell_lines, raw_data and assay_data

## Examples

``` r
get_testdata_codilution()
#> $drug_names
#> [1] "drug_002" "drug_003" "drug_004"
#> 
#> $cell_line_names
#> [1] "cellline_AA" "cellline_BA"
#> 
#> $dt
#>    Drug Name Cell Line Name
#>       <char>         <char>
#> 1:  drug_002    cellline_AA
#> 2:  drug_003    cellline_AA
#> 3:  drug_004    cellline_AA
#> 4:  drug_002    cellline_BA
#> 5:  drug_003    cellline_BA
#> 6:  drug_004    cellline_BA
#>                                                      rId
#>                                                   <char>
#> 1: G00002_drug_002_moa_A_G00001_drug_001_moa_A_72_0.0005
#> 2: G00003_drug_003_moa_A_G00001_drug_001_moa_A_72_0.0005
#> 3: G00004_drug_004_moa_A_G00001_drug_001_moa_A_72_0.0005
#> 4: G00002_drug_002_moa_A_G00001_drug_001_moa_A_72_0.0005
#> 5: G00003_drug_003_moa_A_G00001_drug_001_moa_A_72_0.0005
#> 6: G00004_drug_004_moa_A_G00001_drug_001_moa_A_72_0.0005
#>                                cId   EC50  x_AOC   h RV GEC50 GR Max   ec50
#>                             <char>  <num>  <num> <lgcl> <num>  <num> <lgcl>
#> 1: CL00010_cellline_AA_tissue_x_22 0.9720 0.0280     NA   Inf 0.9720     NA
#> 2: CL00010_cellline_AA_tissue_x_22 0.9700 0.0300     NA   Inf 0.9700     NA
#> 3: CL00010_cellline_AA_tissue_x_22 0.9720 0.0280     NA   Inf 0.9720     NA
#> 4: CL00011_cellline_BA_tissue_x_26 0.9885 0.0115     NA   Inf 0.9885     NA
#> 5: CL00011_cellline_BA_tissue_x_26 0.9935 0.0065     NA   Inf 0.9935     NA
#> 6: CL00011_cellline_BA_tissue_x_26 0.9985 0.0015     NA   Inf 0.9985     NA
#>    GR Inf   GR 0   h GR  E Inf p_value    rss Concentration N_conc    E0
#>    <lgcl> <lgcl> <lgcl> <lgcl>  <lgcl> <lgcl>         <num>  <int> <num>
#> 1:     NA     NA     NA     NA      NA     NA      -3.30103      1     0
#> 2:     NA     NA     NA     NA      NA     NA      -3.30103      1     0
#> 3:     NA     NA     NA     NA      NA     NA      -3.30103      1     0
#> 4:     NA     NA     NA     NA      NA     NA      -3.30103      1     0
#> 5:     NA     NA     NA     NA      NA     NA      -3.30103      1     0
#> 6:     NA     NA     NA     NA      NA     NA      -3.30103      1     0
#>                fit_type normalization_type fit_source Gnumber Drug MOA
#>                  <char>             <char>     <char>  <char>   <char>
#> 1: DRCTooFewPointsToFit                 RV        gDR  G00002    moa_A
#> 2: DRCTooFewPointsToFit                 RV        gDR  G00003    moa_A
#> 3: DRCTooFewPointsToFit                 RV        gDR  G00004    moa_A
#> 4: DRCTooFewPointsToFit                 RV        gDR  G00002    moa_A
#> 5: DRCTooFewPointsToFit                 RV        gDR  G00003    moa_A
#> 6: DRCTooFewPointsToFit                 RV        gDR  G00004    moa_A
#>    Gnumber_2 Drug Name 2 Drug MOA2 Duration Concentration 2    clid   Tissue
#>       <char>      <char>    <char>    <num>           <num>  <char>   <char>
#> 1:    G00001    drug_001     moa_A       72           5e-04 CL00010 tissue_x
#> 2:    G00001    drug_001     moa_A       72           5e-04 CL00010 tissue_x
#> 3:    G00001    drug_001     moa_A       72           5e-04 CL00010 tissue_x
#> 4:    G00001    drug_001     moa_A       72           5e-04 CL00011 tissue_x
#> 5:    G00001    drug_001     moa_A       72           5e-04 CL00011 tissue_x
#> 6:    G00001    drug_001     moa_A       72           5e-04 CL00011 tissue_x
#>    ReferenceDivisionTime   GR50   IC50  E Max GR value
#>                    <num>  <num>  <num>  <num>    <num>
#> 1:                    22 0.9720 0.9720 0.9720   0.9720
#> 2:                    22 0.9700 0.9700 0.9700   0.9700
#> 3:                    22 0.9720 0.9720 0.9720   0.9720
#> 4:                    26 0.9885 0.9885 0.9885   0.9885
#> 5:                    26 0.9935 0.9935 0.9935   0.9935
#> 6:                    26 0.9985 0.9985 0.9985   0.9985
#> 
#> $raw_data
#>                                                                     rId
#>                                                                  <char>
#>   1:              G00002_drug_002_moa_A_G00001_drug_001_moa_A_72_0.0005
#>   2:              G00002_drug_002_moa_A_G00001_drug_001_moa_A_72_0.0005
#>   3: G00002_drug_002_moa_A_G00001_drug_001_moa_A_72_0.00158113883008419
#>   4: G00002_drug_002_moa_A_G00001_drug_001_moa_A_72_0.00158113883008419
#>   5:               G00002_drug_002_moa_A_G00001_drug_001_moa_A_72_0.005
#>  ---                                                                   
#> 104:                 G00004_drug_004_moa_A_G00001_drug_001_moa_A_72_0.5
#> 105:    G00004_drug_004_moa_A_G00001_drug_001_moa_A_72_1.58113883008419
#> 106:    G00004_drug_004_moa_A_G00001_drug_001_moa_A_72_1.58113883008419
#> 107:                   G00004_drug_004_moa_A_G00001_drug_001_moa_A_72_5
#> 108:                   G00004_drug_004_moa_A_G00001_drug_001_moa_A_72_5
#>                                  cId x_mean  x_AOC x_AOC_range  xc50  x_max
#>                               <char>  <num>  <num>      <lgcl> <num>  <num>
#>   1: CL00010_cellline_AA_tissue_x_22 0.9720 0.0280          NA   Inf 0.9720
#>   2: CL00010_cellline_AA_tissue_x_22 0.9827 0.0173          NA   Inf 0.9827
#>   3: CL00010_cellline_AA_tissue_x_22 0.8084 0.1916          NA   Inf 0.8084
#>   4: CL00010_cellline_AA_tissue_x_22 0.8742 0.1258          NA   Inf 0.8742
#>   5: CL00010_cellline_AA_tissue_x_22 0.5138 0.4862          NA   Inf 0.5138
#>  ---                                                                       
#> 104: CL00011_cellline_BA_tissue_x_26 0.5872 0.4128          NA   Inf 0.5872
#> 105: CL00011_cellline_BA_tissue_x_26 0.5272 0.4728          NA   Inf 0.5272
#> 106: CL00011_cellline_BA_tissue_x_26 0.5872 0.4128          NA   Inf 0.5872
#> 107: CL00011_cellline_BA_tissue_x_26 0.5272 0.4728          NA   Inf 0.5272
#> 108: CL00011_cellline_BA_tissue_x_26 0.5872 0.4128          NA   Inf 0.5872
#>        ec50  x_inf    x_0      h     r2 p_value    rss maxlog10Concentration
#>      <lgcl> <lgcl> <lgcl> <lgcl> <lgcl>  <lgcl> <lgcl>                 <num>
#>   1:     NA     NA     NA     NA     NA      NA     NA              -3.30103
#>   2:     NA     NA     NA     NA     NA      NA     NA              -3.30103
#>   3:     NA     NA     NA     NA     NA      NA     NA              -2.80103
#>   4:     NA     NA     NA     NA     NA      NA     NA              -2.80103
#>   5:     NA     NA     NA     NA     NA      NA     NA              -2.30103
#>  ---                                                                        
#> 104:     NA     NA     NA     NA     NA      NA     NA              -0.30103
#> 105:     NA     NA     NA     NA     NA      NA     NA               0.19897
#> 106:     NA     NA     NA     NA     NA      NA     NA               0.19897
#> 107:     NA     NA     NA     NA     NA      NA     NA               0.69897
#> 108:     NA     NA     NA     NA     NA      NA     NA               0.69897
#>      N_conc x_sd_avg             fit_type normalization_type fit_source Gnumber
#>       <int>    <num>               <char>             <char>     <char>  <char>
#>   1:      1        0 DRCTooFewPointsToFit                 RV        gDR  G00002
#>   2:      1        0 DRCTooFewPointsToFit                 GR        gDR  G00002
#>   3:      1        0 DRCTooFewPointsToFit                 RV        gDR  G00002
#>   4:      1        0 DRCTooFewPointsToFit                 GR        gDR  G00002
#>   5:      1        0 DRCTooFewPointsToFit                 RV        gDR  G00002
#>  ---                                                                           
#> 104:      1        0 DRCTooFewPointsToFit                 GR        gDR  G00004
#> 105:      1        0 DRCTooFewPointsToFit                 RV        gDR  G00004
#> 106:      1        0 DRCTooFewPointsToFit                 GR        gDR  G00004
#> 107:      1        0 DRCTooFewPointsToFit                 RV        gDR  G00004
#> 108:      1        0 DRCTooFewPointsToFit                 GR        gDR  G00004
#>      DrugName drug_moa Gnumber_2 DrugName_2 drug_moa_2 Duration Concentration_2
#>        <char>   <char>    <char>     <char>     <char>    <num>           <num>
#>   1: drug_002    moa_A    G00001   drug_001      moa_A       72     0.000500000
#>   2: drug_002    moa_A    G00001   drug_001      moa_A       72     0.000500000
#>   3: drug_002    moa_A    G00001   drug_001      moa_A       72     0.001581139
#>   4: drug_002    moa_A    G00001   drug_001      moa_A       72     0.001581139
#>   5: drug_002    moa_A    G00001   drug_001      moa_A       72     0.005000000
#>  ---                                                                           
#> 104: drug_004    moa_A    G00001   drug_001      moa_A       72     0.500000000
#> 105: drug_004    moa_A    G00001   drug_001      moa_A       72     1.581138830
#> 106: drug_004    moa_A    G00001   drug_001      moa_A       72     1.581138830
#> 107: drug_004    moa_A    G00001   drug_001      moa_A       72     5.000000000
#> 108: drug_004    moa_A    G00001   drug_001      moa_A       72     5.000000000
#>         clid CellLineName   Tissue ReferenceDivisionTime
#>       <char>       <char>   <char>                 <num>
#>   1: CL00010  cellline_AA tissue_x                    22
#>   2: CL00010  cellline_AA tissue_x                    22
#>   3: CL00010  cellline_AA tissue_x                    22
#>   4: CL00010  cellline_AA tissue_x                    22
#>   5: CL00010  cellline_AA tissue_x                    22
#>  ---                                                    
#> 104: CL00011  cellline_BA tissue_x                    26
#> 105: CL00011  cellline_BA tissue_x                    26
#> 106: CL00011  cellline_BA tissue_x                    26
#> 107: CL00011  cellline_BA tissue_x                    26
#> 108: CL00011  cellline_BA tissue_x                    26
#> 
```
