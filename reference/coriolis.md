# Coriolis parameter over a raster

\\f = 2\Omega\sin(\phi)\\, which vanishes at the equator. A grid
spanning it would divide by something arbitrarily close to zero, so the
band where that happens is returned as `NA` rather than as a very large
Rossby number that looks like a real signal.

## Usage

``` r
coriolis(rast)
```

## Arguments

- rast:

  a `SpatRaster`, used for its geometry only

## Value

a `SpatRaster` of `f`, in s^-1
