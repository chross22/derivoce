# Vertical temperature gradient

The difference between surface and bottom temperature in each cell: a
stratification index, large where a warm surface layer sits over cold
deep water and near zero where the column is well mixed.

## Usage

``` r
vertical_gradient(
  env_dat,
  surface = "SST",
  bottom = "BOTT",
  depth = NULL,
  name = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- surface:

  name of the shallower temperature column

- bottom:

  name of the deeper temperature column

- depth:

  optional depth column name; when given, the difference is divided by
  depth to give a per-metre rate rather than a total difference.
  [`datamatch::attach_bathymetry()`](https://camilleross.org/datamatch/reference/attach_bathymetry.html)
  supplies a `DEPTH` column.

- name:

  name for the new column

## Value

`env_dat` with the vertical gradient column added

## Details

Both defaults come from the same Copernicus physics dataset (`SST` and
`BOTT` in the `datamatch` catalog, `thetao` and `bottomT` as Copernicus
codes), so the usual case needs no extra fetch — it is a difference
between two columns already present.

## Other levels, and a better measure

The two columns are arguments, not fixed, so this is not restricted to
the surface and the sea floor. `datamatch::accessCopernicus()` takes a
`depth` argument and returns one level per call, so a second fetch at a
chosen depth gives a temperature at that level, and the difference
between any two levels can be taken the same way.

Where the two levels are known,
[`buoyancy_frequency()`](https://camilleross.org/derivoce/reference/buoyancy_frequency.md)
is the better quantity. This is a temperature difference, which stands
in for stratification only where salinity is uniform — and in the Gulf
of Maine it is not, since Scotian Shelf inflow is fresh enough to
stratify water barely warmer at the surface. \\N^2\\ counts both
contributions with the weights the equation of state gives them.

## See also

[`buoyancy_frequency()`](https://camilleross.org/derivoce/reference/buoyancy_frequency.md),
[`potential_density()`](https://camilleross.org/derivoce/reference/potential_density.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- datamatch::accessCopernicus(vars = c("SST", "BOTT"), ...)
env <- vertical_gradient(env)                 # SST and BOTT are the defaults

# As a per-metre rate, with depth from datamatch's bathymetry
bathy <- datamatch::fetch_bathymetry(bounding_box = bb)
env <- datamatch::attach_bathymetry(env, bathy, "DEPTH")
env <- vertical_gradient(env, depth = "DEPTH")

# Between two chosen levels rather than surface and sea floor
deep <- datamatch::accessCopernicus(vars = "SST", depth = c(90, 100), ...)
env$SST_90m <- deep$SST
env <- vertical_gradient(env, surface = "SST", bottom = "SST_90m")
} # }
```
