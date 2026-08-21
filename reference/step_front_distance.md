# Distance to the nearest frontal cell within one time step

Distance to the nearest frontal cell within one time step

## Usage

``` r
step_front_distance(points, gradient_column, cutoff, scale)
```

## Arguments

- points:

  an `sf` POINT object for a single time step, carrying a gradient
  column

- gradient_column:

  name of the gradient column

- cutoff:

  gradient value at or above which a cell is frontal

- scale:

  divisor converting metres to the requested unit

## Value

numeric vector of distances, one per point
