# Label rotation-dominated patches within one time step

Label rotation-dominated patches within one time step

## Usage

``` r
eddy_patches(rast, u, v, threshold, min_cells, per = "km")
```

## Arguments

- rast:

  a `SpatRaster` carrying the velocity components

- u, v:

  names of the velocity layers

- threshold:

  multiple of sd(W) below which a cell is rotational

- min_cells:

  smallest patch kept

- per:

  `"km"` or `"m"`, for the radius

## Value

a list of `SpatRaster`s: `patches`, `polarity`, `radius`
