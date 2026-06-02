# Remove batch substring from drug id

Gnumber, i.e. "G12345678" is currently the default format of drug_id.
It's also used as a drug name in some cases.

## Usage

``` r
remove_drug_batch(
  drug_vec,
  drug_p = "^G[0-9]{8}",
  sep_p = "[^0-9|^_]",
  batch_p = ".+"
)
```

## Arguments

- drug_vec:

  atomic vector (e.g., character or integer) with drug id(s)

- drug_p:

  string with regex pattern for drug id. Set to Gnumber format by
  default: "G\[0-9\]{8}".

- sep_p:

  string with regex pattern for separator. Set to any character except
  for digit and space

- batch_p:

  string with regex pattern for batch substring. By default set to any
  character(s): ".+"

## Value

charvec with Gnumber(s)

## Details

By default, Gnumber(s) followed by any character (except for underscore
and any digit) and any batch substring are cleaned:

- G00060245.18 =\> G00060245

- G00060245.1-8 =\> G00060245

- G02948263.1-1.DMA =\> G02948263

- Gnumber followed by the codrug

  - G03252046.1-2;G00376771 =\> G03252046

- Gnumber followed by the two codrugs

  - G03256376.1-2;G00376771.1-19;G02557755 =\> G03256376

- Gnumber followed by the drug name

  - G00018838, Cisplatin =\> G00018838

By default, Gnumber(s) followed by the "\_" or digit (regardless the
batch substring) are not cleaned:

- Gnumber with suffix added to prevent duplicated ids

  - G00060245\_(G00060245.1-8)

- too long Gnumber

  - G123456789.1-12

## Examples

``` r
remove_drug_batch("G00060245.18")
#> [1] "G00060245"
remove_drug_batch("G00060245.1-8")
#> [1] "G00060245"
remove_drug_batch("G00060245.1-1.DMA")
#> [1] "G00060245"

remove_drug_batch("G03252046.1-2;G00376771")
#> [1] "G03252046"
remove_drug_batch("G00018838, Cisplatin")
#> [1] "G00018838"
remove_drug_batch("G03256376.1-2;G00376771.1-19;G02557755")
#> [1] "G03256376"
remove_drug_batch("G00060245_(G00060245.1-8)")
#> [1] "G00060245_(G00060245.1-8)"
remove_drug_batch(c("G00060245.18", "G00060245.1-8", "G00060245.1-1.DMA"))
#> [1] "G00060245" "G00060245" "G00060245"

remove_drug_batch("DRUG_01.123", drug_p = "DRUG_[0-9]+")
#> [1] "DRUG_01"
remove_drug_batch("G00001234:22-1", sep_p = ":")
#> [1] "G00001234"
remove_drug_batch("G00001234.28", batch_p = "[0-9]+")
#> [1] "G00001234"
```
