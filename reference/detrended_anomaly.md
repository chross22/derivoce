# Anomaly with the long-term trend removed as well

Uses the same joint fit as
[`decompose_covariate()`](https://camilleross.org/derivoce/reference/decompose_covariate.md)
rather than removing a trend and a climatology one after the other,
because in sequence each absorbs part of the other.

## Usage

``` r
detrended_anomaly(env_dat, values, reference)
```

## Arguments

- env_dat:

  an `sf` POINT object

- values:

  the covariate

- reference:

  `"climatology"` to drop the seasonal cycle too, `"record"` to keep it

## Value

a numeric vector the same length as `values`
