# A stable key for each point's location

Rounded so that coordinates written and re-read at different precisions
still match across time steps.

## Usage

``` r
location_key(env_dat, digits = 6)
```

## Arguments

- env_dat:

  an `sf` POINT object

- digits:

  rounding applied before forming the key

## Value

character vector, one per row
