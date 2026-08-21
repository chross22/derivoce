# Pull a per-time-step index back out as a series

Several functions here compute **one value per time step** and broadcast
it onto every row, so the object keeps its shape and can carry on down a
pipe:
[`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md),
[`section_transport()`](https://camilleross.org/derivoce/reference/section_transport.md),
[`eastern_gom_salinity()`](https://camilleross.org/derivoce/reference/eastern_gom_salinity.md),
[`northeast_channel_inflow()`](https://camilleross.org/derivoce/reference/northeast_channel_inflow.md)
and the rest of the region-scale indices all behave this way. That is
right for modelling, where the covariate has to line up with the
observations, and wrong for almost everything else. Plotting a 22-year
monthly index from a broadcast column means plotting each value a few
thousand times; writing one out means exporting a file mostly made of
repetition.

## Usage

``` r
index_series(env_dat, vars = NULL)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- vars:

  index columns to extract, or `NULL` to take every column that is
  constant within each time step

## Value

a data frame with `YEAR`, `MONTH`, `DAY` and one column per index, one
row per time step, ordered by date. Not an `sf` object: an index has no
location

## Details

This collapses them back: one row per time step, in date order, as a
plain data frame.

## What it will not do

A column that varies within a time step is a map, not an index, and
collapsing it would silently throw away the spatial pattern and keep an
arbitrary one of its values. Naming such a column is an error rather
than a quiet mean. If a summary of a map per step is what you want, that
is a different operation and an explicit one — take the mean yourself,
or use
[`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md),
which is exactly that with a region attached.

With `vars = NULL` the constant-within-step columns are found for you,
so passing an object through several index functions and then calling
this returns all of them at once.

## See also

[`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md),
[`section_transport()`](https://camilleross.org/derivoce/reference/section_transport.md),
[`derived_indices()`](https://camilleross.org/derivoce/reference/derived_indices.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- eastern_gom_salinity(env)
env <- section_transport(env, from = c(-67.5, 44.5), to = c(-66.0, 43.5))

series <- index_series(env)
plot(with(series, as.Date(paste(YEAR, MONTH, DAY, sep = "-"))),
     series$egom_salinity, type = "l")
} # }
```
