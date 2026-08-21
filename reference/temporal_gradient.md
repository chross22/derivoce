# Temporal gradient of a covariate

Rate of change between consecutive time steps at each location — how
fast conditions are shifting, as distinct from what they currently are.
A water mass warming rapidly is a different habitat from one sitting at
the same temperature.

## Usage

``` r
temporal_gradient(
  env_dat,
  vars = NULL,
  per = c("step", "day", "month"),
  suffix = "_tgrad"
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- vars:

  covariate columns to differentiate; `NULL` does all of them

- per:

  time unit for the rate: `"step"` (default, change per time step),
  `"day"`, or `"month"`

- suffix:

  suffix for the new columns

## Value

`env_dat` with a `<var>_tgrad` column per covariate

## Details

The first time step has no predecessor and is `NA`.

## Examples

``` r
if (FALSE) { # \dontrun{
env <- temporal_gradient(env, "SST")                 # change per time step
env <- temporal_gradient(env, "SST", per = "day")    # degrees C per day
} # }
```
