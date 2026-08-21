# Apply a raster-valued function to every time step

Handles the split/rasterize/compute/join cycle that each spatial
derivative shares, so those functions only have to say what to do to one
raster.

## Usage

``` r
per_time_step(env_dat, vars, fun)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- vars:

  covariate columns the function needs

- fun:

  a function taking a `SpatRaster` and returning a `SpatRaster` whose
  layer names become the new columns

## Value

`env_dat` with the returned layers added as columns
