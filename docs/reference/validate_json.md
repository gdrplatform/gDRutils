# Validate JSON against a schema.

Validate JSON describing an object against a schema.

## Usage

``` r
validate_json(json, schema_path)
```

## Arguments

- json:

  String of JSON in memory.

- schema_path:

  String of the schema to validate against.

## Value

Boolean of whether or not JSON successfully validated.

## Details

This is most often used to validate JSON before passing it in as a
document to an ElasticSearch index.

## Examples

``` r
json <- '{}'
```
