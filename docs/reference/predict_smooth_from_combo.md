# Predict a smoothed response for a drug combination

Predicts a 'smooth' value for a single pair of drug concentrations by
snapping to the nearest available models from a metrics table and
averaging their predictions. This is the combination equivalent of
'predict_efficacy_from_conc'.

## Usage

``` r
predict_smooth_from_combo(conc_1, conc_2, metrics_merged)
```

## Arguments

- conc_1:

  A single numeric value for the desired concentration of the first
  drug.

- conc_2:

  A single numeric value for the desired concentration of the second
  drug.

- metrics_merged:

  A data.table containing all pre-calculated curve fit parameters.
  Expects columns: 'dilution_drug', 'cotrt_value', 'ratio', 'ec50', 'h',
  'x_inf', 'x_0'.

## Value

A single numeric value for the predicted 'smooth' response.

## Examples

``` r
mae <- get_synthetic_data("combo_matrix")
se <- mae[[gDRutils::get_supported_experiments("combo")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se[1, 1], "Metrics")[normalization_type == "RV"]
predict_smooth_from_combo(conc_1 = 1.2, conc_2 = 9.8, metrics_merged = dt_metrics)
#> Requested: (1.20, 9.80) ==> Using models for nearest concentrations: (1.00, 10.00)
#> [1] 0.3168861
```
