# Distance to isobaths

Convenience wrapper on
[`distance_to_contour()`](https://camilleross.org/derivoce/reference/distance_to_contour.md)
for depth contours, the most common use.

## Usage

``` r
distance_to_isobath(
  env_dat,
  depth = "DEPTH",
  levels = c(50, 100, 200),
  per = c("km", "m")
)
```

## Arguments

- env_dat:

  an `sf` POINT object carrying a depth column.
  [`datamatch::attach_bathymetry()`](https://camilleross.org/datamatch/reference/attach_bathymetry.html)
  adds one, named `DEPTH`.

- depth:

  name of the depth column, as a positive magnitude

- levels:

  isobaths to measure to, in the units of `depth`

- per:

  distance unit for the result

## Value

`env_dat` with one distance column per isobath

## Examples

``` r
if (FALSE) { # \dontrun{
bathy <- datamatch::fetch_bathymetry(bounding_box = bb)
env <- datamatch::attach_bathymetry(env, bathy, "DEPTH")
env <- distance_to_isobath(env, levels = c(50, 100, 200))
} # }
```
