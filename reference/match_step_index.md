# Index of each row's time step within the ordered step table

Index of each row's time step within the ordered step table

## Usage

``` r
match_step_index(env_dat, steps)
```

## Arguments

- env_dat:

  an `sf` POINT object

- steps:

  a time-step table from
  [`time_steps()`](https://camilleross.org/derivoce/reference/time_steps.md)

## Value

integer vector, one per row of `env_dat`
