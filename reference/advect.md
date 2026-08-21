# Advect particles through a time-varying velocity field

Fourth-order Runge-Kutta in lon/lat, with velocities converted from m/s
to degrees per day at each evaluation point.

## Usage

``` r
advect(positions, start_time, duration, step_days, velocity, times)
```

## Arguments

- positions:

  two-column matrix of starting lon/lat

- start_time:

  start time, in days

- duration:

  signed integration length in days; negative runs backwards

- step_days:

  step size in days

- velocity:

  list of two-layer `SpatRaster`s, one per time step

- times:

  numeric time of each raster, in days

## Value

two-column matrix of final lon/lat; `NA` for particles that left the
domain or ran aground
