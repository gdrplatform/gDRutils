# Process and save a batch of results.

Process and save a batch of results.

## Usage

``` r
process_batch(
  batch,
  start_index,
  fun_name,
  unique_id,
  total_iterations,
  temp_dir,
  FUN,
  ...
)
```

## Arguments

- batch:

  A subset of the vector or list `x` to be processed.

- start_index:

  Integer indicating the starting index of the batch in the original
  vector `x`.

- fun_name:

  Character string representing the name of the function `FUN` for use
  in file naming.

- unique_id:

  String with unique identifier for the current task and user to ensure
  file uniqueness.

- total_iterations:

  Integer indicating the total number of iterations in the original
  vector `x`.

- temp_dir:

  Character string specifying the directory where batch results are
  saved.

- FUN:

  A user-defined function to apply to each element of the batch.

- ...:

  Optional arguments passed to `FUN`.

## Value

This function does not return a value. It saves the processed batch
results to disk as a `.qs2` file.

## Details

The function applies `FUN` to each element in `batch`, saves the results
to a file named according to the format
`<fun_name>_<unique_id>_<start_index>_of_<total_iterations>_batch.qs2`,
and clears memory using [`gc()`](https://rdrr.io/r/base/gc.html) after
saving.

## Examples

``` r
process_batch(list(1, 2, 3), 100, "my_function", "unique_task_id_user", 1000, tempdir(), function(x) x^2)
```
