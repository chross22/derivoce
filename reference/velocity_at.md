# Velocity in degrees per day at arbitrary positions and time

Bilinear in space, linear in time between the two bracketing fields.
Sampling the two bracketing rasters at the particle positions and
blending the values is equivalent to blending the rasters first, and
much cheaper.

## Usage

``` r
velocity_at(positions, time, velocity, times)
```

## Arguments

- positions:

  two-column matrix of lon/lat

- time:

  time in days

- velocity:

  list of two-layer `SpatRaster`s

- times:

  numeric time of each raster, in days

## Value

two-column matrix of dlon/dt and dlat/dt, in degrees per day
