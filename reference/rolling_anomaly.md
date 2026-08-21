# Anomaly against a centred rolling mean

Anomaly against a centred rolling mean

## Usage

``` r
rolling_anomaly(env_dat, values, location, window)
```

## Arguments

- env_dat:

  an `sf` POINT object

- values:

  the covariate values

- location:

  per-row location key

- window:

  number of time steps in the rolling mean

## Value

numeric vector of anomalies
