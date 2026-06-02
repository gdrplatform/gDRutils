# Update gDR synonyms for the identifier

Update gDR synonyms for the identifier

## Usage

``` r
update_idfs_synonyms(data, dict = get_idfs_synonyms())
```

## Arguments

- data:

  list of charvec with identifiers data

- dict:

  list with dictionary

## Value

list

## Examples

``` r
mdict <- list(duration = "time")
iv <- c("Time", "Duration", "time")
update_idfs_synonyms(iv, dict = mdict)
#> [1] "Duration" "Duration" "Duration"
```
