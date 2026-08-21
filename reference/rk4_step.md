# One fourth-order Runge-Kutta step

One fourth-order Runge-Kutta step

## Usage

``` r
rk4_step(positions, time, dt, velocity, times)
```

## Arguments

- positions:

  two-column matrix of lon/lat

- time:

  current time in days

- dt:

  signed step size in days

- velocity:

  list of two-layer `SpatRaster`s

- times:

  numeric time of each raster, in days

## Value

two-column matrix of updated lon/lat
