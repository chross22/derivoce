# Cells a contour passes through

A contour almost never falls exactly on a cell centre, so a cell is
marked when the level lies between its own value and any of its
neighbours' — that is, when the field crosses the level somewhere in the
cell's neighbourhood.

## Usage

``` r
contour_cells(grid, level)
```

## Arguments

- grid:

  a one-layer `SpatRaster`

- level:

  the value being contoured

## Value

a `SpatRaster` of `1` on the contour and `NA` elsewhere, as
[`terra::distance()`](https://rspatial.github.io/terra/reference/distance.html)
expects
