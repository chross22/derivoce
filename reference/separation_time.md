# Time for a parcel pair to reach a target separation

Advects each seed together with an eastward and a northward companion,
checking after every step whether either pair has separated far enough.

## Usage

``` r
separation_time(
  seeds,
  delta_0,
  delta_f,
  start_time,
  sign,
  max_days,
  step_days,
  velocity,
  times
)
```

## Arguments

- seeds:

  two-column matrix of lon/lat

- delta_0:

  initial separation, in km

- delta_f:

  target separation, in km

- start_time:

  start time in days

- sign:

  `-1` to integrate backwards, `1` forwards

- max_days:

  give-up time

- step_days:

  step size in days

- velocity:

  list of two-layer `SpatRaster`s

- times:

  numeric time of each raster, in days

## Value

numeric vector of times in days; `NA` where the target was never met
