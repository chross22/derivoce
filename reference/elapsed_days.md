# Days elapsed since the first time step

Real dates rather than step positions, so a gap in the record leaves a
gap in the trend rather than compressing it.

## Usage

``` r
elapsed_days(env_dat)
```

## Arguments

- env_dat:

  an `sf` POINT object

## Value

numeric vector of days, one per row
