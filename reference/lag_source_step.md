# Which earlier step each step draws its lagged value from

Returns one index per time step, or `NA` where the lag lands outside the
record. Separating this from
[`lag_covariate()`](https://camilleross.org/derivoce/reference/lag_covariate.md)
keeps the calendar arithmetic in one place and lets it be tested on the
step table alone.

## Usage

``` r
lag_source_step(steps, n, by)
```

## Arguments

- steps:

  a time-step table from
  [`time_steps()`](https://camilleross.org/derivoce/reference/time_steps.md)

- n:

  how far back, in `by` units

- by:

  `"step"`, `"day"`, `"month"`, or `"year"`

## Value

integer vector of source step indices, `NA` where there is none
