# Split environmental data into its time steps

Split environmental data into its time steps

## Usage

``` r
time_steps(env_dat)
```

## Arguments

- env_dat:

  an `sf` POINT object with `YEAR`/`MONTH`/`DAY` columns

## Value

a data frame of unique time steps, ordered chronologically
