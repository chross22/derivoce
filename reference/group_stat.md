# Apply a summary within each group and broadcast it back

[`stats::ave()`](https://rdrr.io/r/stats/ave.html) would drop rows where
the covariate is `NA`, so the summary is taken explicitly with `na.rm`
and matched back by group.

## Usage

``` r
group_stat(values, groups, fun)
```

## Arguments

- values:

  a numeric vector

- groups:

  a grouping vector of the same length

- fun:

  a summary function

## Value

a numeric vector the same length as `values`
