# Check whether coordinates are evenly spaced

The tolerance has to clear the quantisation in the coordinates
themselves. Copernicus stores longitude and latitude as **float32**,
whose spacing near 67 degrees resolves to about 8e-6, so a nominally
uniform 1/12-degree grid arrives with spacings varying by around 1e-4 of
a cell. That is a property of the file format, not of the grid, and a
tolerance tight enough to reject it rejects every real Copernicus
download.

## Usage

``` r
is_regular(values, tolerance = 0.001)
```

## Arguments

- values:

  sorted unique coordinate values

- tolerance:

  relative tolerance on the spacing

## Value

`TRUE` if the spacing is constant within tolerance

## Details

1e-3 is loose enough for that and still far tighter than any genuine
irregularity. Scattered observations have spacings that differ by order
one, and a grid missing a row or column has one interval of double
width, so both are still caught by a wide margin.

The comparison is against the median spacing rather than the first, so a
single odd interval at the start cannot set the reference for everything
else.
