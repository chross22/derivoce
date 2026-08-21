# FTLE from a flow map

Differentiates the flow map to get the deformation gradient, forms the
Cauchy-Green strain tensor, and returns the exponential separation rate.

## Usage

``` r
ftle_from_flow_map(seeds, final, integration_days)
```

## Arguments

- seeds:

  two-column matrix of starting lon/lat, on a regular grid

- final:

  two-column matrix of ending lon/lat

- integration_days:

  integration length in days, unsigned

## Value

numeric vector of FTLE values, in 1/day

## Details

Positions are converted to metres on an equirectangular plane about the
grid centre before differentiating, so the strain tensor is computed in
a single consistent metric rather than mixing degrees of longitude and
latitude, which have different physical lengths.
