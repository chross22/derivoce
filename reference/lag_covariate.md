# Lagged covariate values

Adds each covariate's value from `n` units earlier at the same location.
Populations respond to conditions with a delay. A bloom feeds the
animals sampled a month later, not the ones sampled during it, so the
lagged value can carry more signal than the concurrent one.

## Usage

``` r
lag_covariate(
  env_dat,
  vars = NULL,
  n = 1,
  by = c("step", "day", "month", "year"),
  suffix = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- vars:

  covariate columns to lag; `NULL` does all of them

- n:

  how far to look back, counted in `by` units. A vector adds one column
  per lag.

- by:

  what `n` counts: `"step"` (the default), `"day"`, `"month"`, or
  `"year"`

- suffix:

  suffix for the new columns. The default includes `n`, and the unit too
  unless it is `"step"`, so `SST_lag1` and `SST_lag1month` cannot be
  confused for each other. Supply one entry per lag when `n` has
  several.

## Value

`env_dat` with a lagged column per covariate and lag. Steps with no
predecessor are `NA`.

## Lag by calendar units, not by position

`by` decides what `n` counts, and the distinction matters:

- `"step"` (the default) counts **positions in the series**. `n = 1` is
  the previous time step, whatever period that represents. This is only
  unambiguous when the steps are evenly spaced and complete.

- `"day"`, `"month"`, `"year"` count **calendar time**. `n = 3` with
  `by = "month"` finds the step stamped exactly three calendar months
  earlier, and returns `NA` if there is no such step.

Prefer a calendar unit whenever the lag has a biological meaning. "Three
months ago" is a statement about the organism; "three steps ago" is a
statement about how the data happened to be fetched, and the two stop
agreeing the moment a month is missing from the record.

A gap makes them disagree silently. In a monthly series missing April,
`by = "step"` treats March as May's predecessor, so a one-step lag
quietly becomes a two-month one. `by = "month"` returns `NA` for May
instead, because April genuinely is not there.

Calendar lags are matched on the exact `YEAR`/`MONTH`/`DAY` stamp. For
monthly products, whose day is always 1, that is exact. For daily
products, `by = "month"` from the 31st looks for a 31st, and months that
have no 31st return `NA` rather than silently sliding to the 30th.

Locations are matched by coordinate, so this assumes a fixed grid across
time steps, which is what gridded products give.

## Reproducing the published lag

Ross et al. (2023) used a **one-month** lag of sea surface temperature.
On a complete monthly series that is `n = 1` either way, but the
faithful form is the calendar one:

    lag_covariate(env, "SST", n = 1, by = "month")

The two part company the moment a month is missing from the record,
where `by = "step"` reaches back to whatever step precedes the gap and
calls it one month. Use `by = "month"` when the intent is the published
lag.

## Several lags at once

`n` may be a vector, which adds one column per lag. This is what an
autoregressive design needs: a covariate's own past at a series of
offsets, entered together as predictors.

    lag_covariate(env, "SST", n = 1:3, by = "year")
    # adds SST_lag1year, SST_lag2year, SST_lag3year

With `by = "year"` this gives the same calendar month in each preceding
year, so the seasonal cycle is held fixed and what remains is the
interannual signal. That is usually the intended comparison, and it is
not what `by = "step"` with `n = 12` gives on a record with any month
missing.

## References

Ross C, Runge J, Roberts J, Brady D, Tupper B, Record N (2023).
Estimating North Atlantic right whale prey based on Calanus finmarchicus
thresholds. *Marine Ecology Progress Series* **703**, 1-16.
[doi:10.3354/meps14204](https://doi.org/10.3354/meps14204)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- lag_covariate(env, "SST")                      # SST_lag1, previous step
env <- lag_covariate(env, "CHL", n = 2)               # CHL_lag2

# Calendar lags, which say what they mean
env <- lag_covariate(env, "CHL", n = 3, by = "month") # CHL_lag3month
env <- lag_covariate(env, "SST", n = 1, by = "year")  # same month last year
env <- lag_covariate(env, "SST", n = 30, by = "day")  # daily products

# An autoregressive set: the same month in each of the last three years
env <- lag_covariate(env, "SST", n = 1:3, by = "year")
} # }
```
