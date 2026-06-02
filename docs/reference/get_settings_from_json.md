# Get settings from JSON file In most common scenario the settings are stored in JSON file to avoid hardcoding

Get settings from JSON file In most common scenario the settings are
stored in JSON file to avoid hardcoding

## Usage

``` r
get_settings_from_json(
  s = NULL,
  json_path = system.file(package = "gDRutils", "settings.json")
)
```

## Arguments

- s:

  charvec with setting entry/entries

- json_path:

  string with the path to the JSON file

## Value

value/values for entry/entries from JSON file

## Examples

``` r
if (!nchar(system.file(package="gDRutils"))) {
   get_settings_from_json()
}
```
