# Departure of a covariate from its mean, per location

Departure of a covariate from its mean, per location

## Usage

``` r
anomaly(env_dat, var, reference)
```

## Arguments

- env_dat:

  an `sf` POINT object

- var:

  covariate column name

- reference:

  `"record"`, `"climatology"`, or a rolling window length

## Value

numeric vector of anomalies, one per row
