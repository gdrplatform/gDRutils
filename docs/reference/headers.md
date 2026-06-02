# Get or reset headers for one or all header field(s) respectively

Get the expected header(s) for one field or reset all header fields

## Usage

``` r
get_header(k = NULL)
```

## Arguments

- k:

  string of field (data type) to return headers for

## Value

For `get_header` a character vector of headers for field `k`.

## Details

If `get_header` is called with no values, the entire available header
list is returned.

## Examples

``` r
get_header(k = NULL)
#> $manifest
#> $manifest$barcode
#> [1] "Barcode" "Plate"  
#> 
#> $manifest$template
#> [1] "Template"  "Treatment"
#> 
#> $manifest$duration
#> [1] "Duration"
#> 
#> 
#> $raw_data
#> [1] "ReadoutValue"    "BackgroundValue" "UntrtReadout"    "Day0Readout"    
#> [5] "masked"         
#> 
#> $normalized_results
#> [1] "x"                    "CorrectedReadout"     "GRvalue"             
#> [4] "RelativeViability"    "DivisionTime"         "RefGRvalue"          
#> [7] "RefRelativeViability"
#> 
#> $averaged_results
#> [1] "x"                     "x_std"                 "std_RelativeViability"
#> [4] "std_GRvalue"           "count"                 "x_sd"                 
#> [7] "x_std_sd"             
#> 
#> $response_metrics
#>  [1] "x_mean"                "x_AOC"                 "x_AOC_range"          
#>  [4] "xc50"                  "x_max"                 "ec50"                 
#>  [7] "x_inf"                 "x_0"                   "h"                    
#> [10] "r2"                    "p_value"               "rss"                  
#> [13] "maxlog10Concentration" "N_conc"                "x_sd_avg"             
#> [16] "fit_type"              "x_mean_sd"             "x_AOC_sd"             
#> [19] "x_AOC_range_sd"        "xc50_sd"               "x_max_sd"             
#> [22] "ec50_sd"               "x_inf_sd"              "x_0_sd"               
#> [25] "h_sd"                  "r2_sd"                 "x_sd_avg_sd"          
#> 
#> $metric_average_fields
#> $metric_average_fields$mean
#> [1] "x_mean"      "x_AOC"       "x_AOC_range" "x_max"       "x_inf"      
#> [6] "x_0"        
#> 
#> $metric_average_fields$geometric_mean
#>  [1] "xc50"    "ec50"    "GR50"    "GEC50"   "IC50"    "EC50"    "GR_xc50"
#>  [8] "RV_xc50" "GR_ec50" "RV_ec50"
#> 
#> $metric_average_fields$fit_type
#> [1] "fit_type"    "Fit Type"    "Fit Type RV" "Fit Type GR" "RV_fit_type"
#> [6] "GR_fit_type"
#> 
#> $metric_average_fields$blacklisted
#> [1] "cellline_tissue"              "Tissue"                      
#> [3] "cellline_ref_div_time"        "Reference Division Time"     
#> [5] "ReferenceDivisionTime"        "cellline_parental_identifier"
#> [7] "Parental Identifier"          "parental_identifier"         
#> 
#> 
#> $add_clid
#> $add_clid$cellline_name
#> [1] "CellLineName"
#> 
#> $add_clid$cellline_tissue
#> [1] "Tissue"
#> 
#> $add_clid$cellline_parental_identifier
#> [1] "parental_identifier"
#> 
#> $add_clid$cellline_subtype
#> [1] "subtype"
#> 
#> $add_clid$cellline_ref_div_time
#> [1] "ReferenceDivisionTime"
#> 
#> 
#> $metrics_names
#>    x_mean    x_AOC    x_AOC_range    xc50   x_max    ec50    x_inf    x_0   
#> RV "RV_mean" "RV_AOC" "RV_AOC_range" "IC50" "E_max"  "EC50"  "E_inf"  "E_0" 
#> GR "GR_mean" "GR_AOC" "GR_AOC_range" "GR50" "GR_max" "GEC50" "GR_inf" "GR_0"
#>    h      r2      p_value      rss      maxlog10Concentration      N_conc     
#> RV "h_RV" "RV_r2" "RV_p_value" "RV_rss" "RV_maxlog10Concentration" "RV_N_conc"
#> GR "h_GR" "GR_r2" "GR_p_value" "GR_rss" "GR_maxlog10Concentration" "GR_N_conc"
#>    x_sd_avg    fit_type      x_mean_sd    x_AOC_sd    x_AOC_range_sd   
#> RV "RV_sd_avg" "fit_type_RV" "RV_mean_sd" "RV_AOC_sd" "RV_AOC_range_sd"
#> GR "GR_sd_avg" "fit_type_GR" "GR_mean_sd" "GR_AOC_sd" "GR_AOC_range_sd"
#>    xc50_sd   x_max_sd    ec50_sd    x_inf_sd    x_0_sd    h_sd      r2_sd     
#> RV "IC50_sd" "E_max_sd"  "EC50_sd"  "E_inf_sd"  "E_0_sd"  "h_RV_sd" "RV_r2_sd"
#> GR "GR50_sd" "GR_max_sd" "GEC50_sd" "GR_inf_sd" "GR_0_sd" "h_GR_sd" "GR_r2_sd"
#>    x_sd_avg_sd   
#> RV "RV_sd_avg_sd"
#> GR "GR_sd_avg_sd"
#> 
#> $metrics_results
#>  [1] "maxlog10Concentration"    "maxlog10Concentration_sd"
#>  [3] "N_conc"                   "N_conc_sd"               
#>  [5] "cotrt_value"              "cotrt_value_sd"          
#>  [7] "ratio"                    "ratio_sd"                
#>  [9] "dilution_drug"            "count"                   
#> [11] "x_mean"                   "x_AOC"                   
#> [13] "x_AOC_range"              "xc50"                    
#> [15] "x_max"                    "ec50"                    
#> [17] "x_inf"                    "x_0"                     
#> [19] "h"                        "r2"                      
#> [21] "p_value"                  "rss"                     
#> [23] "maxlog10Concentration"    "N_conc"                  
#> [25] "x_sd_avg"                 "fit_type"                
#> [27] "x_mean_sd"                "x_AOC_sd"                
#> [29] "x_AOC_range_sd"           "xc50_sd"                 
#> [31] "x_max_sd"                 "ec50_sd"                 
#> [33] "x_inf_sd"                 "x_0_sd"                  
#> [35] "h_sd"                     "r2_sd"                   
#> [37] "x_sd_avg_sd"              "RV_mean"                 
#> [39] "GR_mean"                  "RV_AOC"                  
#> [41] "GR_AOC"                   "RV_AOC_range"            
#> [43] "GR_AOC_range"             "IC50"                    
#> [45] "GR50"                     "E_max"                   
#> [47] "GR_max"                   "EC50"                    
#> [49] "GEC50"                    "E_inf"                   
#> [51] "GR_inf"                   "E_0"                     
#> [53] "GR_0"                     "h_RV"                    
#> [55] "h_GR"                     "RV_r2"                   
#> [57] "GR_r2"                    "RV_p_value"              
#> [59] "GR_p_value"               "RV_rss"                  
#> [61] "GR_rss"                   "RV_maxlog10Concentration"
#> [63] "GR_maxlog10Concentration" "RV_N_conc"               
#> [65] "GR_N_conc"                "RV_sd_avg"               
#> [67] "GR_sd_avg"                "fit_type_RV"             
#> [69] "fit_type_GR"              "RV_mean_sd"              
#> [71] "GR_mean_sd"               "RV_AOC_sd"               
#> [73] "GR_AOC_sd"                "RV_AOC_range_sd"         
#> [75] "GR_AOC_range_sd"          "IC50_sd"                 
#> [77] "GR50_sd"                  "E_max_sd"                
#> [79] "GR_max_sd"                "EC50_sd"                 
#> [81] "GEC50_sd"                 "E_inf_sd"                
#> [83] "GR_inf_sd"                "E_0_sd"                  
#> [85] "GR_0_sd"                  "h_RV_sd"                 
#> [87] "h_GR_sd"                  "RV_r2_sd"                
#> [89] "GR_r2_sd"                 "RV_sd_avg_sd"            
#> [91] "GR_sd_avg_sd"            
#> 
#> $controlled
#> $controlled[[1]]
#> [1] "clid"
#> 
#> $controlled$barcode
#> [1] "Barcode" "Plate"  
#> 
#> $controlled$template
#> [1] "Template"  "Treatment"
#> 
#> $controlled$duration
#> [1] "Duration"
#> 
#> $controlled[[5]]
#> [1] "Gnumber"
#> 
#> $controlled[[6]]
#> [1] "Concentration"
#> 
#> $controlled[[7]]
#> [1] "Gnumber_2"
#> 
#> $controlled[[8]]
#> [1] "Gnumber_3"
#> 
#> $controlled[[9]]
#> [1] "Gnumber_4"
#> 
#> $controlled[[10]]
#> [1] "Gnumber_5"
#> 
#> $controlled[[11]]
#> [1] "Gnumber_6"
#> 
#> $controlled[[12]]
#> [1] "Gnumber_7"
#> 
#> $controlled[[13]]
#> [1] "Gnumber_8"
#> 
#> $controlled[[14]]
#> [1] "Gnumber_9"
#> 
#> $controlled[[15]]
#> [1] "Gnumber_10"
#> 
#> $controlled[[16]]
#> [1] "Concentration_2"
#> 
#> $controlled[[17]]
#> [1] "Concentration_3"
#> 
#> $controlled[[18]]
#> [1] "Concentration_4"
#> 
#> $controlled[[19]]
#> [1] "Concentration_5"
#> 
#> $controlled[[20]]
#> [1] "Concentration_6"
#> 
#> $controlled[[21]]
#> [1] "Concentration_7"
#> 
#> $controlled[[22]]
#> [1] "Concentration_8"
#> 
#> $controlled[[23]]
#> [1] "Concentration_9"
#> 
#> $controlled[[24]]
#> [1] "Concentration_10"
#> 
#> 
#> $reserved
#> $reserved$cellline_name
#> [1] "CellLineName"
#> 
#> $reserved$cellline_tissue
#> [1] "Tissue"
#> 
#> $reserved$cellline_parental_identifier
#> [1] "parental_identifier"
#> 
#> $reserved$cellline_subtype
#> [1] "subtype"
#> 
#> $reserved$cellline_ref_div_time
#> [1] "ReferenceDivisionTime"
#> 
#> $reserved[[6]]
#> [1] "DrugName"
#> 
#> $reserved[[7]]
#> [1] "masked"
#> 
#> $reserved[[8]]
#> [1] "DrugName_2"
#> 
#> $reserved[[9]]
#> [1] "DrugName_3"
#> 
#> $reserved[[10]]
#> [1] "DrugName_4"
#> 
#> $reserved[[11]]
#> [1] "DrugName_5"
#> 
#> $reserved[[12]]
#> [1] "DrugName_6"
#> 
#> $reserved[[13]]
#> [1] "DrugName_7"
#> 
#> $reserved[[14]]
#> [1] "DrugName_8"
#> 
#> $reserved[[15]]
#> [1] "DrugName_9"
#> 
#> $reserved[[16]]
#> [1] "DrugName_10"
#> 
#> $reserved[[17]]
#> [1] "drug_moa_2"
#> 
#> $reserved[[18]]
#> [1] "drug_moa_3"
#> 
#> $reserved[[19]]
#> [1] "drug_moa_4"
#> 
#> $reserved[[20]]
#> [1] "drug_moa_5"
#> 
#> $reserved[[21]]
#> [1] "drug_moa_6"
#> 
#> $reserved[[22]]
#> [1] "drug_moa_7"
#> 
#> $reserved[[23]]
#> [1] "drug_moa_8"
#> 
#> $reserved[[24]]
#> [1] "drug_moa_9"
#> 
#> $reserved[[25]]
#> [1] "drug_moa_10"
#> 
#> $reserved[[26]]
#> [1] "ReadoutValue"
#> 
#> $reserved[[27]]
#> [1] "BackgroundValue"
#> 
#> $reserved[[28]]
#> [1] "UntrtReadout"
#> 
#> $reserved[[29]]
#> [1] "Day0Readout"
#> 
#> $reserved[[30]]
#> [1] "masked"
#> 
#> $reserved[[31]]
#> [1] "x"
#> 
#> $reserved[[32]]
#> [1] "CorrectedReadout"
#> 
#> $reserved[[33]]
#> [1] "GRvalue"
#> 
#> $reserved[[34]]
#> [1] "RelativeViability"
#> 
#> $reserved[[35]]
#> [1] "DivisionTime"
#> 
#> $reserved[[36]]
#> [1] "RefGRvalue"
#> 
#> $reserved[[37]]
#> [1] "RefRelativeViability"
#> 
#> $reserved[[38]]
#> [1] "x"
#> 
#> $reserved[[39]]
#> [1] "x_std"
#> 
#> $reserved[[40]]
#> [1] "std_RelativeViability"
#> 
#> $reserved[[41]]
#> [1] "std_GRvalue"
#> 
#> $reserved[[42]]
#> [1] "count"
#> 
#> $reserved[[43]]
#> [1] "x_sd"
#> 
#> $reserved[[44]]
#> [1] "x_std_sd"
#> 
#> $reserved[[45]]
#> [1] "maxlog10Concentration"
#> 
#> $reserved[[46]]
#> [1] "maxlog10Concentration_sd"
#> 
#> $reserved[[47]]
#> [1] "N_conc"
#> 
#> $reserved[[48]]
#> [1] "N_conc_sd"
#> 
#> $reserved[[49]]
#> [1] "cotrt_value"
#> 
#> $reserved[[50]]
#> [1] "cotrt_value_sd"
#> 
#> $reserved[[51]]
#> [1] "ratio"
#> 
#> $reserved[[52]]
#> [1] "ratio_sd"
#> 
#> $reserved[[53]]
#> [1] "dilution_drug"
#> 
#> $reserved[[54]]
#> [1] "count"
#> 
#> $reserved[[55]]
#> [1] "x_mean"
#> 
#> $reserved[[56]]
#> [1] "x_AOC"
#> 
#> $reserved[[57]]
#> [1] "x_AOC_range"
#> 
#> $reserved[[58]]
#> [1] "xc50"
#> 
#> $reserved[[59]]
#> [1] "x_max"
#> 
#> $reserved[[60]]
#> [1] "ec50"
#> 
#> $reserved[[61]]
#> [1] "x_inf"
#> 
#> $reserved[[62]]
#> [1] "x_0"
#> 
#> $reserved[[63]]
#> [1] "h"
#> 
#> $reserved[[64]]
#> [1] "r2"
#> 
#> $reserved[[65]]
#> [1] "p_value"
#> 
#> $reserved[[66]]
#> [1] "rss"
#> 
#> $reserved[[67]]
#> [1] "maxlog10Concentration"
#> 
#> $reserved[[68]]
#> [1] "N_conc"
#> 
#> $reserved[[69]]
#> [1] "x_sd_avg"
#> 
#> $reserved[[70]]
#> [1] "fit_type"
#> 
#> $reserved[[71]]
#> [1] "x_mean_sd"
#> 
#> $reserved[[72]]
#> [1] "x_AOC_sd"
#> 
#> $reserved[[73]]
#> [1] "x_AOC_range_sd"
#> 
#> $reserved[[74]]
#> [1] "xc50_sd"
#> 
#> $reserved[[75]]
#> [1] "x_max_sd"
#> 
#> $reserved[[76]]
#> [1] "ec50_sd"
#> 
#> $reserved[[77]]
#> [1] "x_inf_sd"
#> 
#> $reserved[[78]]
#> [1] "x_0_sd"
#> 
#> $reserved[[79]]
#> [1] "h_sd"
#> 
#> $reserved[[80]]
#> [1] "r2_sd"
#> 
#> $reserved[[81]]
#> [1] "x_sd_avg_sd"
#> 
#> $reserved[[82]]
#> [1] "RV_mean"
#> 
#> $reserved[[83]]
#> [1] "GR_mean"
#> 
#> $reserved[[84]]
#> [1] "RV_AOC"
#> 
#> $reserved[[85]]
#> [1] "GR_AOC"
#> 
#> $reserved[[86]]
#> [1] "RV_AOC_range"
#> 
#> $reserved[[87]]
#> [1] "GR_AOC_range"
#> 
#> $reserved[[88]]
#> [1] "IC50"
#> 
#> $reserved[[89]]
#> [1] "GR50"
#> 
#> $reserved[[90]]
#> [1] "E_max"
#> 
#> $reserved[[91]]
#> [1] "GR_max"
#> 
#> $reserved[[92]]
#> [1] "EC50"
#> 
#> $reserved[[93]]
#> [1] "GEC50"
#> 
#> $reserved[[94]]
#> [1] "E_inf"
#> 
#> $reserved[[95]]
#> [1] "GR_inf"
#> 
#> $reserved[[96]]
#> [1] "E_0"
#> 
#> $reserved[[97]]
#> [1] "GR_0"
#> 
#> $reserved[[98]]
#> [1] "h_RV"
#> 
#> $reserved[[99]]
#> [1] "h_GR"
#> 
#> $reserved[[100]]
#> [1] "RV_r2"
#> 
#> $reserved[[101]]
#> [1] "GR_r2"
#> 
#> $reserved[[102]]
#> [1] "RV_p_value"
#> 
#> $reserved[[103]]
#> [1] "GR_p_value"
#> 
#> $reserved[[104]]
#> [1] "RV_rss"
#> 
#> $reserved[[105]]
#> [1] "GR_rss"
#> 
#> $reserved[[106]]
#> [1] "RV_maxlog10Concentration"
#> 
#> $reserved[[107]]
#> [1] "GR_maxlog10Concentration"
#> 
#> $reserved[[108]]
#> [1] "RV_N_conc"
#> 
#> $reserved[[109]]
#> [1] "GR_N_conc"
#> 
#> $reserved[[110]]
#> [1] "RV_sd_avg"
#> 
#> $reserved[[111]]
#> [1] "GR_sd_avg"
#> 
#> $reserved[[112]]
#> [1] "fit_type_RV"
#> 
#> $reserved[[113]]
#> [1] "fit_type_GR"
#> 
#> $reserved[[114]]
#> [1] "RV_mean_sd"
#> 
#> $reserved[[115]]
#> [1] "GR_mean_sd"
#> 
#> $reserved[[116]]
#> [1] "RV_AOC_sd"
#> 
#> $reserved[[117]]
#> [1] "GR_AOC_sd"
#> 
#> $reserved[[118]]
#> [1] "RV_AOC_range_sd"
#> 
#> $reserved[[119]]
#> [1] "GR_AOC_range_sd"
#> 
#> $reserved[[120]]
#> [1] "IC50_sd"
#> 
#> $reserved[[121]]
#> [1] "GR50_sd"
#> 
#> $reserved[[122]]
#> [1] "E_max_sd"
#> 
#> $reserved[[123]]
#> [1] "GR_max_sd"
#> 
#> $reserved[[124]]
#> [1] "EC50_sd"
#> 
#> $reserved[[125]]
#> [1] "GEC50_sd"
#> 
#> $reserved[[126]]
#> [1] "E_inf_sd"
#> 
#> $reserved[[127]]
#> [1] "GR_inf_sd"
#> 
#> $reserved[[128]]
#> [1] "E_0_sd"
#> 
#> $reserved[[129]]
#> [1] "GR_0_sd"
#> 
#> $reserved[[130]]
#> [1] "h_RV_sd"
#> 
#> $reserved[[131]]
#> [1] "h_GR_sd"
#> 
#> $reserved[[132]]
#> [1] "RV_r2_sd"
#> 
#> $reserved[[133]]
#> [1] "GR_r2_sd"
#> 
#> $reserved[[134]]
#> [1] "RV_sd_avg_sd"
#> 
#> $reserved[[135]]
#> [1] "GR_sd_avg_sd"
#> 
#> $reserved[[136]]
#> [1] "WellRow"
#> 
#> $reserved[[137]]
#> [1] "WellColumn"
#> 
#> 
#> $ordered_1
#>  [1] "CellLineName"     "Tissue"           "Duration"         "DrugName"        
#>  [5] "Concentration"    "DrugName_2"       "Concentration_2"  "DrugName_3"      
#>  [9] "Concentration_3"  "DrugName_4"       "Concentration_4"  "DrugName_5"      
#> [13] "Concentration_5"  "DrugName_6"       "Concentration_6"  "DrugName_7"      
#> [17] "Concentration_7"  "DrugName_8"       "Concentration_8"  "DrugName_9"      
#> [21] "Concentration_9"  "DrugName_10"      "Concentration_10"
#> 
#> $ordered_2
#> $ordered_2[[1]]
#> [1] "x"
#> 
#> $ordered_2[[2]]
#> [1] "CorrectedReadout"
#> 
#> $ordered_2[[3]]
#> [1] "GRvalue"
#> 
#> $ordered_2[[4]]
#> [1] "RelativeViability"
#> 
#> $ordered_2[[5]]
#> [1] "DivisionTime"
#> 
#> $ordered_2[[6]]
#> [1] "RefGRvalue"
#> 
#> $ordered_2[[7]]
#> [1] "RefRelativeViability"
#> 
#> $ordered_2[[8]]
#> [1] "x"
#> 
#> $ordered_2[[9]]
#> [1] "x_std"
#> 
#> $ordered_2[[10]]
#> [1] "std_RelativeViability"
#> 
#> $ordered_2[[11]]
#> [1] "std_GRvalue"
#> 
#> $ordered_2[[12]]
#> [1] "count"
#> 
#> $ordered_2[[13]]
#> [1] "x_sd"
#> 
#> $ordered_2[[14]]
#> [1] "x_std_sd"
#> 
#> $ordered_2[[15]]
#> [1] "maxlog10Concentration"
#> 
#> $ordered_2[[16]]
#> [1] "maxlog10Concentration_sd"
#> 
#> $ordered_2[[17]]
#> [1] "N_conc"
#> 
#> $ordered_2[[18]]
#> [1] "N_conc_sd"
#> 
#> $ordered_2[[19]]
#> [1] "cotrt_value"
#> 
#> $ordered_2[[20]]
#> [1] "cotrt_value_sd"
#> 
#> $ordered_2[[21]]
#> [1] "ratio"
#> 
#> $ordered_2[[22]]
#> [1] "ratio_sd"
#> 
#> $ordered_2[[23]]
#> [1] "dilution_drug"
#> 
#> $ordered_2[[24]]
#> [1] "count"
#> 
#> $ordered_2[[25]]
#> [1] "x_mean"
#> 
#> $ordered_2[[26]]
#> [1] "x_AOC"
#> 
#> $ordered_2[[27]]
#> [1] "x_AOC_range"
#> 
#> $ordered_2[[28]]
#> [1] "xc50"
#> 
#> $ordered_2[[29]]
#> [1] "x_max"
#> 
#> $ordered_2[[30]]
#> [1] "ec50"
#> 
#> $ordered_2[[31]]
#> [1] "x_inf"
#> 
#> $ordered_2[[32]]
#> [1] "x_0"
#> 
#> $ordered_2[[33]]
#> [1] "h"
#> 
#> $ordered_2[[34]]
#> [1] "r2"
#> 
#> $ordered_2[[35]]
#> [1] "p_value"
#> 
#> $ordered_2[[36]]
#> [1] "rss"
#> 
#> $ordered_2[[37]]
#> [1] "maxlog10Concentration"
#> 
#> $ordered_2[[38]]
#> [1] "N_conc"
#> 
#> $ordered_2[[39]]
#> [1] "x_sd_avg"
#> 
#> $ordered_2[[40]]
#> [1] "fit_type"
#> 
#> $ordered_2[[41]]
#> [1] "x_mean_sd"
#> 
#> $ordered_2[[42]]
#> [1] "x_AOC_sd"
#> 
#> $ordered_2[[43]]
#> [1] "x_AOC_range_sd"
#> 
#> $ordered_2[[44]]
#> [1] "xc50_sd"
#> 
#> $ordered_2[[45]]
#> [1] "x_max_sd"
#> 
#> $ordered_2[[46]]
#> [1] "ec50_sd"
#> 
#> $ordered_2[[47]]
#> [1] "x_inf_sd"
#> 
#> $ordered_2[[48]]
#> [1] "x_0_sd"
#> 
#> $ordered_2[[49]]
#> [1] "h_sd"
#> 
#> $ordered_2[[50]]
#> [1] "r2_sd"
#> 
#> $ordered_2[[51]]
#> [1] "x_sd_avg_sd"
#> 
#> $ordered_2[[52]]
#> [1] "RV_mean"
#> 
#> $ordered_2[[53]]
#> [1] "GR_mean"
#> 
#> $ordered_2[[54]]
#> [1] "RV_AOC"
#> 
#> $ordered_2[[55]]
#> [1] "GR_AOC"
#> 
#> $ordered_2[[56]]
#> [1] "RV_AOC_range"
#> 
#> $ordered_2[[57]]
#> [1] "GR_AOC_range"
#> 
#> $ordered_2[[58]]
#> [1] "IC50"
#> 
#> $ordered_2[[59]]
#> [1] "GR50"
#> 
#> $ordered_2[[60]]
#> [1] "E_max"
#> 
#> $ordered_2[[61]]
#> [1] "GR_max"
#> 
#> $ordered_2[[62]]
#> [1] "EC50"
#> 
#> $ordered_2[[63]]
#> [1] "GEC50"
#> 
#> $ordered_2[[64]]
#> [1] "E_inf"
#> 
#> $ordered_2[[65]]
#> [1] "GR_inf"
#> 
#> $ordered_2[[66]]
#> [1] "E_0"
#> 
#> $ordered_2[[67]]
#> [1] "GR_0"
#> 
#> $ordered_2[[68]]
#> [1] "h_RV"
#> 
#> $ordered_2[[69]]
#> [1] "h_GR"
#> 
#> $ordered_2[[70]]
#> [1] "RV_r2"
#> 
#> $ordered_2[[71]]
#> [1] "GR_r2"
#> 
#> $ordered_2[[72]]
#> [1] "RV_p_value"
#> 
#> $ordered_2[[73]]
#> [1] "GR_p_value"
#> 
#> $ordered_2[[74]]
#> [1] "RV_rss"
#> 
#> $ordered_2[[75]]
#> [1] "GR_rss"
#> 
#> $ordered_2[[76]]
#> [1] "RV_maxlog10Concentration"
#> 
#> $ordered_2[[77]]
#> [1] "GR_maxlog10Concentration"
#> 
#> $ordered_2[[78]]
#> [1] "RV_N_conc"
#> 
#> $ordered_2[[79]]
#> [1] "GR_N_conc"
#> 
#> $ordered_2[[80]]
#> [1] "RV_sd_avg"
#> 
#> $ordered_2[[81]]
#> [1] "GR_sd_avg"
#> 
#> $ordered_2[[82]]
#> [1] "fit_type_RV"
#> 
#> $ordered_2[[83]]
#> [1] "fit_type_GR"
#> 
#> $ordered_2[[84]]
#> [1] "RV_mean_sd"
#> 
#> $ordered_2[[85]]
#> [1] "GR_mean_sd"
#> 
#> $ordered_2[[86]]
#> [1] "RV_AOC_sd"
#> 
#> $ordered_2[[87]]
#> [1] "GR_AOC_sd"
#> 
#> $ordered_2[[88]]
#> [1] "RV_AOC_range_sd"
#> 
#> $ordered_2[[89]]
#> [1] "GR_AOC_range_sd"
#> 
#> $ordered_2[[90]]
#> [1] "IC50_sd"
#> 
#> $ordered_2[[91]]
#> [1] "GR50_sd"
#> 
#> $ordered_2[[92]]
#> [1] "E_max_sd"
#> 
#> $ordered_2[[93]]
#> [1] "GR_max_sd"
#> 
#> $ordered_2[[94]]
#> [1] "EC50_sd"
#> 
#> $ordered_2[[95]]
#> [1] "GEC50_sd"
#> 
#> $ordered_2[[96]]
#> [1] "E_inf_sd"
#> 
#> $ordered_2[[97]]
#> [1] "GR_inf_sd"
#> 
#> $ordered_2[[98]]
#> [1] "E_0_sd"
#> 
#> $ordered_2[[99]]
#> [1] "GR_0_sd"
#> 
#> $ordered_2[[100]]
#> [1] "h_RV_sd"
#> 
#> $ordered_2[[101]]
#> [1] "h_GR_sd"
#> 
#> $ordered_2[[102]]
#> [1] "RV_r2_sd"
#> 
#> $ordered_2[[103]]
#> [1] "GR_r2_sd"
#> 
#> $ordered_2[[104]]
#> [1] "RV_sd_avg_sd"
#> 
#> $ordered_2[[105]]
#> [1] "GR_sd_avg_sd"
#> 
#> $ordered_2[[106]]
#> [1] "ReadoutValue"
#> 
#> $ordered_2[[107]]
#> [1] "BackgroundValue"
#> 
#> $ordered_2[[108]]
#> [1] "UntrtReadout"
#> 
#> $ordered_2[[109]]
#> [1] "Day0Readout"
#> 
#> $ordered_2[[110]]
#> [1] "masked"
#> 
#> $ordered_2[[111]]
#> [1] "ReferenceDivisionTime"
#> 
#> $ordered_2[[112]]
#> [1] "clid"
#> 
#> $ordered_2[[113]]
#> [1] "Gnumber"
#> 
#> $ordered_2[[114]]
#> [1] "Gnumber_2"
#> 
#> $ordered_2[[115]]
#> [1] "Gnumber_3"
#> 
#> $ordered_2[[116]]
#> [1] "Gnumber_4"
#> 
#> $ordered_2[[117]]
#> [1] "Gnumber_5"
#> 
#> $ordered_2[[118]]
#> [1] "Gnumber_6"
#> 
#> $ordered_2[[119]]
#> [1] "Gnumber_7"
#> 
#> $ordered_2[[120]]
#> [1] "Gnumber_8"
#> 
#> $ordered_2[[121]]
#> [1] "Gnumber_9"
#> 
#> $ordered_2[[122]]
#> [1] "Gnumber_10"
#> 
#> $ordered_2$barcode
#> [1] "Barcode" "Plate"  
#> 
#> $ordered_2$template
#> [1] "Template"  "Treatment"
#> 
#> $ordered_2$duration
#> [1] "Duration"
#> 
#> $ordered_2[[126]]
#> [1] "WellRow"
#> 
#> $ordered_2[[127]]
#> [1] "WellColumn"
#> 
#> 
#> $id
#> [1] "rId" "cId"
#> 
#> $iso_position
#> [1] "iso_level" "pos_x"     "pos_y"     "pos_x_ref" "pos_y_ref"
#> 
#> $excess
#> [1] "smooth"       "hsa_excess"   "bliss_excess"
#> 
#> $excess_results
#> [1] "smooth"          "hsa_excess"      "bliss_excess"    "smooth_sd"      
#> [5] "hsa_excess_sd"   "bliss_excess_sd"
#> 
#> $scores
#> [1] "hsa_score"   "bliss_score" "CIScore_50"  "CIScore_80" 
#> 
#> $scores_results
#> [1] "hsa_score"      "bliss_score"    "CIScore_50"     "CIScore_80"    
#> [5] "hsa_score_sd"   "bliss_score_sd" "CIScore_50_sd"  "CIScore_80_sd" 
#> 
#> $isobolograms
#> [1] "normalization_type" "iso_level"          "pos_x"             
#> [4] "pos_y"              "pos_x_ref"          "pos_y_ref"         
#> [7] "log2_CI"            "log10_ratio_conc"  
#> 
#> $isobolograms_results
#>  [1] "normalization_type"    "iso_level"             "pos_x"                
#>  [4] "pos_y"                 "pos_x_ref"             "pos_y_ref"            
#>  [7] "log2_CI"               "log10_ratio_conc"      "normalization_type_sd"
#> [10] "iso_level_sd"          "pos_x_sd"              "pos_y_sd"             
#> [13] "pos_x_ref_sd"          "pos_y_ref_sd"          "log2_CI_sd"           
#> [16] "log10_ratio_conc_sd"  
#> 
#> $fit_source
#> [1] "fit_source"
#> 
#> $obsolete
#> [1] "RV"     "GR"     "Excess"
#> 
get_header("manifest")
#> $barcode
#> [1] "Barcode" "Plate"  
#> 
#> $template
#> [1] "Template"  "Treatment"
#> 
#> $duration
#> [1] "Duration"
#> 
```
