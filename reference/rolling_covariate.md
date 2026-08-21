# Rolling summary of a covariate over a trailing window

Summarises each location's recent history: the mean of the last three
months, the variability of the last year, the coldest step in the last
six. Where
[`integrate_covariate()`](https://camilleross.org/derivoce/reference/integrate_covariate.md)
accumulates a total, this describes the distribution the total came
from, which is often the more useful covariate — a mean and a standard
deviation say different things about a place than their sum does.

## Usage

``` r
rolling_covariate(
  env_dat,
  vars = NULL,
  n = 3,
  by = c("step", "day", "month", "year"),
  stat = c("mean", "sd", "min", "max", "sum", "median", "range"),
  min_obs = 1L,
  suffix = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- vars:

  covariate columns, or `NULL` for all numeric ones

- n:

  length of the window, in `by` units, including the current step

- by:

  `"step"` to count positions in the record, or `"day"`, `"month"`,
  `"year"` to count calendar time

- stat:

  one or more of `"mean"`, `"sd"`, `"min"`, `"max"`, `"sum"`,
  `"median"`, `"range"`

- min_obs:

  fewest non-missing values a window must hold to be summarised

- suffix:

  one per `stat`, or `NULL` to name them automatically

## Value

`env_dat` with one column per covariate per statistic

## Details

The window is **trailing and inclusive**: it ends at the current step
and includes it. A three-month mean at March covers January, February
and March.

## Steps or calendar time

`by = "step"` counts positions in the record and `by = "day"`, `"month"`
or `"year"` count calendar time, exactly as in
[`lag_covariate()`](https://camilleross.org/derivoce/reference/lag_covariate.md).
The two agree until the record has a gap and then disagree silently: on
a monthly series missing April, a three-*step* window at June covers
March, May and June, while a three-*month* window covers April, May and
June and finds only two of them.

Which is right depends on the question. "The mean of the last three
months" is a statement about the ocean and wants `by = "month"`. "The
mean of the last three observations" is a statement about the record.

## Windows that are not full

Early steps have less history behind them than the window asks for, and
a location can be absent from some of the steps in it. `min_obs` sets
how many values a window must actually contain before it is summarised;
below that the result is `NA` rather than a mean of whatever happened to
be there. The default of 1 is permissive, so the first steps of a record
get a summary of a short window rather than nothing. Raise it if a
partial window would mislead.

## See also

[`integrate_covariate()`](https://camilleross.org/derivoce/reference/integrate_covariate.md),
[`lag_covariate()`](https://camilleross.org/derivoce/reference/lag_covariate.md),
[`cell_anomaly()`](https://camilleross.org/derivoce/reference/cell_anomaly.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Conditions over the season leading up to each observation.
env <- rolling_covariate(env, "SST", n = 3, by = "month")

# How variable it has been, which is a different covariate from how warm.
env <- rolling_covariate(env, "SST", n = 12, by = "month", stat = "sd")
} # }
```
