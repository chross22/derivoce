# Which steps fall in the trailing window ending at a given step

The calendar cases are done on the same month counter as
[`lag_source_step()`](https://camilleross.org/derivoce/reference/lag_source_step.md),
so that a window and a lag of the same size agree about what a month is.

## Usage

``` r
window_steps(steps, i, n, by)
```

## Arguments

- steps:

  a time-step table from
  [`time_steps()`](https://camilleross.org/derivoce/reference/time_steps.md)

- i:

  the step the window ends at

- n:

  window length in `by` units, including step `i`

- by:

  `"step"`, `"day"`, `"month"`, or `"year"`

## Value

the indices of the contributing steps
