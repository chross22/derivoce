# Columns that are bookkeeping rather than data

Kept identical to `datamatch::is_bookkeeping_column()`. Duplicated
rather than imported for the same reason
[`covariate_columns()`](https://camilleross.org/derivoce/reference/covariate_columns.md)
is: this package works on the shape datamatch returns, not on datamatch.

## Usage

``` r
is_bookkeeping_column(names)
```

## Arguments

- names:

  column names to inspect

## Value

one per name
