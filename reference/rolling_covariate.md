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
  by,
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

  what `n` counts: `"hour"`, `"day"`, `"month"` or `"year"`. Required.

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
and March. On a record finer than monthly, `"month"` and `"year"`
windows reach back by calendar month but stop at the current step, so a
daily record's monthly maximum on the 10th covers the 1st to the 10th,
never the days after it.

## Calendar time

`by` is required and says what `n` counts: `"hour"`, `"day"`, `"month"`
or `"year"`, as in
[`lag_covariate()`](https://camilleross.org/derivoce/reference/lag_covariate.md).
There is no option to count positions in the record and no default unit.
On a monthly series missing April, a three-*step* window at June would
cover March, May and June and call it three months. A calendar window
covers April, May and June and finds only two of them, which `min_obs`
can then reject.

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
