# Horizontal spatial gradient of a covariate

Adds the magnitude of each covariate's horizontal gradient — how fast it
changes with distance — which is what picks out fronts, the convergence
zones where plankton aggregate.

## Usage

``` r
horizontal_gradient(
  env_dat,
  vars = NULL,
  per = c("km", "m"),
  components = FALSE,
  suffix = "_grad"
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- vars:

  covariate columns to differentiate; `NULL` does all of them

- per:

  distance unit for the result: `"km"` (default) or `"m"`

- components:

  also add the eastward and northward components, as `<var>_grad_x` and
  `<var>_grad_y`

- suffix:

  suffix for the magnitude column

## Value

`env_dat` with a `<var>_grad` column per covariate, in covariate units
per `per`

## Details

Computed by central differences on the covariate's own grid, in real
distance units rather than degrees. That distinction matters: a degree
of longitude is about 83 km at 42 degrees N but 111 km at the equator,
so a gradient left in per-degree units is stretched by latitude and not
comparable across a study area.

This is a deliberate departure from `raster::terrain()`, which the
original pipeline of Ross et al. (2023) used for `sst_grad` and
`uv_grad`. `terrain()` treats its input as an elevation in the same
units as the coordinates and returns a slope angle, which is
dimensionally meaningless for a field measured in degrees Celsius. The
magnitude here is a real rate of change with a real unit.

## References

Ross C, Runge J, Roberts J, Brady D, Tupper B, Record N (2023).
Estimating North Atlantic right whale prey based on Calanus finmarchicus
thresholds. *Marine Ecology Progress Series* **703**, 1-16.
[doi:10.3354/meps14204](https://doi.org/10.3354/meps14204)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- datamatch::accessCopernicus(vars = c("SST", "SSS"), ...)
env <- horizontal_gradient(env, "SST")        # SST_grad, in degrees C per km
env <- horizontal_gradient(env, "SST", components = TRUE)

# `vars = NULL` does every covariate column, which is right for a plain fetch
# and wrong once static or non-numeric columns have been attached.
env <- horizontal_gradient(env)
} # }
```
