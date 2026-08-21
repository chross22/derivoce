# Sample the velocity field at positions and a time

Sample the velocity field at positions and a time

## Usage

``` r
sample_velocity(positions, time, velocity, times)
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

two-column matrix of eastward and northward velocity, in m/s
