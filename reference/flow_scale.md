# Typical flow speed and domain size, for diagnosing empty Lagrangian output

FTLE and FSLE both fail the same way: a domain that is small relative to
how far the flow carries a parcel leaves nothing to report. The numbers
needed to say so are the median speed and the extent, so they are
gathered once here.

## Usage

``` r
flow_scale(env_dat, u, v)
```

## Arguments

- env_dat:

  an `sf` POINT object

- u, v:

  velocity column names, in m/s

## Value

list with `speed` (m/s), `width_km`, and `height_km`
