# Anomaly of a covariate averaged over a box

Averages a covariate over a lon/lat box in each time step, then
subtracts a reference to give an anomaly. The result is **one number per
time step**, broadcast to every row, so it behaves like a climate index
rather than a map.

## Usage

``` r
box_anomaly(
  env_dat,
  var,
  box,
  reference = c("climatology", "record", "none"),
  name = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- var:

  covariate column to average

- box:

  named list with `xmin`, `xmax`, `ymin`, `ymax`, in degrees

- reference:

  `"climatology"` (the default) removes a separate mean per calendar
  month, so only departures from the usual conditions for that month
  survive. `"record"` removes one mean over the whole series, leaving
  the seasonal cycle in. `"none"` returns the box mean itself.

- name:

  name for the new column

## Value

`env_dat` with an anomaly column added, in the units of `var`

## Details

This is the simplest of the region-scale indices and often the most
robust. It asks whether conditions in a place were unusual, without
needing velocities or endmembers, so it survives on products where the
other methods cannot run.

## What it cannot tell you

A box mean says conditions changed, not that water moved. A fresh
anomaly in the eastern Gulf of Maine is consistent with more Scotian
Shelf inflow, and also with local runoff, rainfall, or ice melt. Where a
transport across a section measures the crossing directly, this measures
its most visible consequence and asks you to supply the interpretation.

## See also

[`eastern_gom_salinity()`](https://camilleross.org/derivoce/reference/eastern_gom_salinity.md),
[`derived_indices()`](https://camilleross.org/derivoce/reference/derived_indices.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- box_anomaly(env, "SSS", box = list(
  xmin = -68.5, xmax = -66.5, ymin = 43.0, ymax = 44.5
))
} # }
```
