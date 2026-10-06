# Temporal gradient of a covariate

Rate of change since the previous calendar `by`, at each location: how
fast conditions are shifting, as distinct from what they currently are.
A water mass warming rapidly is a different habitat from one sitting at
the same temperature.

## Usage

``` r
temporal_gradient(env_dat, vars = NULL, by, per = by, suffix = "_tgrad")
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- vars:

  covariate columns to differentiate; `NULL` does all of them

- by:

  the calendar unit to look back by: `"hour"`, `"day"`, `"month"` or
  `"year"`. Required, with no default, for the reasons given in
  [`lag_covariate()`](https://camilleross.org/derivoce/reference/lag_covariate.md).

- per:

  time unit the rate is expressed in: `"hour"`, `"day"`, `"month"` or
  `"year"`. Defaults to `by`, so a monthly lag gives the change per
  month exactly. A different unit divides by the actual elapsed days,
  using a mean month of 30.4375 days and a year of 365.25.

- suffix:

  suffix for the new columns

## Value

`env_dat` with a `<var>_tgrad` column per covariate

## Details

The "previous" step is the one stamped exactly one `by` earlier, as in
[`lag_covariate()`](https://camilleross.org/derivoce/reference/lag_covariate.md),
so a gap gives `NA` and not a rate that quietly spans two months. A step
with no such predecessor, including the first, is `NA`.

## Examples

``` r
if (FALSE) { # \dontrun{
env <- temporal_gradient(env, "SST", by = "month")                 # per month
env <- temporal_gradient(env, "SST", by = "month", per = "day")    # degrees C per day
} # }
```
