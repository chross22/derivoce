# Rasterize one time step's points onto their own grid

Spatial derivatives are only defined on a grid, so the points have to be
put back onto one. Gridded ocean products come as a regular lon/lat
lattice, and datamatch's access functions flatten that lattice to points
without moving them, so the grid can be recovered exactly from the
unique coordinates. Copernicus, HYCOM, CCMP and most ERDDAP grids are of
that kind.

## Usage

``` r
rasterize_step(points, vars)
```

## Arguments

- points:

  an `sf` POINT object for a single time step

- vars:

  covariate columns to include as layers

## Value

a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with one layer per variable

## Details

An unstructured mesh is not.
[`datamatch::accessFVCOM()`](https://camilleross.org/datamatch/reference/accessFVCOM.html)
returns one row per mesh node, and those nodes are irregularly spaced by
design, which is the point of a mesh: resolution follows the coastline
rather than a lattice.

Irregular points are rejected rather than interpolated, whether they
come from a mesh or from scattered observations. Silently gridding them
would produce a gradient field that looks plausible and is mostly
interpolation artifact. The derivations that need no lattice – every
temporal one, plus
[`distance_to_shore()`](https://camilleross.org/derivoce/reference/distance_to_shore.md)
and
[`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md)
– work on a mesh unchanged.
