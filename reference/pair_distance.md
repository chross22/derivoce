# Distance between paired positions, in km

Equirectangular about a reference latitude, which is accurate for the
small separations FSLE deals with and consistent with the metric used
elsewhere here.

## Usage

``` r
pair_distance(a, b, reference_lat)
```

## Arguments

- a, b:

  two-column matrices of lon/lat

- reference_lat:

  latitude to evaluate the longitude scaling at

## Value

numeric vector of distances in km
