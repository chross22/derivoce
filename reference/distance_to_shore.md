# Distance to the nearest shoreline

Distance from each point to the nearest coast, in kilometres. Static —
it does not vary in time — so it is computed once for the unique
locations and shared across every time step.

## Usage

``` r
distance_to_shore(
  env_dat,
  resolution = c("medium", "large", "small"),
  margin = 5,
  name = "shore_dist"
)
```

## Source

Coastlines are Natural Earth, via `rnaturalearth`
(<https://www.naturalearthdata.com/>). Public domain, and the
`resolution` argument selects between their 1:110m, 1:50m, and 1:10m
physical coastline layers.

The resolution is a real choice rather than a detail. A deeply indented
coast like Maine's is shortened by a coarse layer, so distances near
shore come out too large. `"large"` is the 1:10m layer and needs
`rnaturalearthhires`.

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- resolution:

  Natural Earth coastline resolution: `"medium"`, `"large"`, or
  `"small"`

- margin:

  degrees of padding around the data when cropping the coastline. Land
  just outside the study area can still be the nearest shore, so
  cropping tightly to the bounding box would overstate distances at the
  edges.

- name:

  name for the new column

## Value

`env_dat` with a distance-to-shore column, in km

## Details

Distance from shore is a broad proxy for several things at once: depth,
terrestrial nutrient input, tidal mixing, and larval retention all
covary with it. That makes it a useful covariate and a poor explanation
— a model leaning on it is telling you *where*, not *why*. This is the
covariate the older pipeline had as `dist`.

## Choosing a coastline

`resolution` selects the Natural Earth coastline detail:

- `"medium"` (the default, 1:50m) suits shelf-scale work.

- `"large"` (1:10m) resolves individual islands and inlets, which
  matters in a place like the Gulf of Maine where the coast is deeply
  indented, but is slower and needs `rnaturalearthhires`.

- `"small"` (1:110m) is too coarse for coastal work; it smooths bays
  away entirely.

Distances are computed against true geodesic distance on the ellipsoid
rather than a projected approximation, so they stay correct across a
wide domain.

## Examples

``` r
if (FALSE) { # \dontrun{
env <- distance_to_shore(env)
# Finer coastline, for a deeply indented shore:
env <- distance_to_shore(env, resolution = "large")
} # }
```
