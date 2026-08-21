# Distance to a contour of a covariate

Distance from every point to the nearest place where a covariate crosses
a given value. Position relative to a contour is often what matters
rather than the value itself — plankton distributions track the shelf
break, and "20 km inshore of the 100 m isobath" locates that far better
than "depth = 85 m" does.

## Usage

``` r
distance_to_contour(env_dat, var, levels, per = c("km", "m"), prefix = NULL)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- var:

  covariate whose contours are wanted

- levels:

  one or more values to contour

- per:

  distance unit for the result: `"km"` (default) or `"m"`

- prefix:

  stem for the new column names; defaults to `var`

## Value

`env_dat` with one distance column per level. A level lying entirely
outside the covariate's range in a time step gives `NA` for that step.

## Details

A cell is treated as on the contour if `level` falls between its value
and that of any neighbour, so the contour is found even though no cell
sits exactly on it.

## Examples

``` r
if (FALSE) { # \dontrun{
# Distance to the shelf break and two other isobaths
env <- distance_to_contour(env, "DEPTH", levels = c(50, 100, 200))
# -> DEPTH_dist_50, DEPTH_dist_100, DEPTH_dist_200

# Works on any field: distance to the 10 degree isotherm
env <- distance_to_contour(env, "SST", levels = 10)
} # }
```
