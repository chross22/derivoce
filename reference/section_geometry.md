# Sample points, unit normal, and segment length for a section

Geometry is done on an equirectangular plane about the section's mean
latitude, the same metric the Lagrangian code uses, so a degree of
longitude and a degree of latitude are comparable lengths before any
normal is taken.

## Usage

``` r
section_geometry(env_dat, from, to, spacing = NULL)
```

## Arguments

- env_dat:

  an `sf` POINT object, used for the grid spacing

- from, to:

  endpoints, `c(longitude, latitude)`

- spacing:

  sample spacing in km, or `NULL` for half a grid cell

## Value

list with `points` (matrix of lon/lat), `normal` (unit vector), and `ds`
(segment length in metres)
