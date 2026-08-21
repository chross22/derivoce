# Warn when there is not enough history to form the requested anomaly

Warn when there is not enough history to form the requested anomaly

## Usage

``` r
warn_thin_climatology(groups, reference, standardize)
```

## Arguments

- groups:

  the grouping vector

- reference:

  `"climatology"` or `"record"`

- standardize:

  whether a standard deviation is also needed

## Value

invisible `NULL`, called for the warning
