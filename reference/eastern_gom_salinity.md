# Eastern Gulf of Maine salinity index

Surface salinity anomaly over the eastern Gulf of Maine, the region
where Scotian Shelf inflow first appears after rounding Cape Sable.
Negative values are fresher than usual, which indicates a stronger
inflow of cold fresh Scotian Shelf Water.

## Usage

``` r
eastern_gom_salinity(
  env_dat,
  salinity = "SSS",
  reference = c("climatology", "record", "none"),
  name = "egom_salinity"
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- salinity:

  name of the surface salinity column, in PSU

- reference:

  passed to
  [`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md)

- name:

  name for the new column

## Value

`env_dat` with an `egom_salinity` column, in PSU

## Details

Follows the approach of Grodsky et al. (2025), who showed that satellite
surface salinity in the eastern Gulf tracks winter Scotian Shelf inflow
and relates it to the coastal and interior pathways the water then
takes. Their index is built from SMAP satellite salinity; this computes
the same kind of quantity from whatever gridded salinity field you
supply, so it is **not a reproduction of their published series** and
should not be compared to it value for value.

## Why this is the seasonal one

The signal is a winter one. Scotian Shelf Water arrives in the upper
layers of the Gulf in winter, roughly six months after the late-summer
deep influx of warm salty water. A whole-year mean mixes the two regimes
together and dilutes the thing being measured, so restrict to winter
months at the fetch, or subset afterwards, rather than reading an annual
value as an inflow index.

## References

Grodsky SA, Vandemark D, Levin J (2025). An eastern Gulf of Maine
salinity index for monitoring winter Scotian Shelf inflow and its
relation to coastal and interior pathways. *Journal of Geophysical
Research: Oceans* **130**(5).
[doi:10.1029/2024JC021891](https://doi.org/10.1029/2024JC021891)

## See also

[`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md),
[`eastern_gom_box()`](https://camilleross.org/derivoce/reference/eastern_gom_box.md),
[`scotian_shelf_inflow()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- eastern_gom_salinity(env)
} # }
```
