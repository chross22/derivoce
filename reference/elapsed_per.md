# Elapsed time between each step and the one a calendar unit earlier

Elapsed time between each step and the one a calendar unit earlier

## Usage

``` r
elapsed_per(steps, by, per)
```

## Arguments

- steps:

  a time-step table from
  [`time_steps()`](https://camilleross.org/derivoce/reference/time_steps.md)

- by:

  the calendar unit one step looks back by

- per:

  the unit the elapsed time is expressed in

## Value

numeric vector, one per step; `NA` where there is no such step. It is
exactly 1 when `per` equals `by`.
