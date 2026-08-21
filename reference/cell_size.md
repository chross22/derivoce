# Cell size in real distance units

Converts a lon/lat grid's degree spacing into distance. Longitude
spacing shrinks with the cosine of latitude, so it is computed per row
rather than once for the whole grid — at 45 degrees the error from
ignoring that is about 30 percent.

## Usage

``` r
cell_size(layer, per = "km")
```

## Arguments

- layer:

  a `SpatRaster`

- per:

  `"km"` or `"m"`

## Value

`list(x =, y =)`, where `x` is a `SpatRaster` of per-row longitude
spacing and `y` a scalar

## Details

A projected raster is already in linear units, so its resolution is used
directly.
