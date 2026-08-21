# Buoyancy frequency between two depths

How strongly the water column is stratified: the frequency at which a
parcel displaced vertically would oscillate, squared. Written \\N^2\\,
and equal to \\-(g/\rho_0)\\\partial\rho/\partial z\\.

## Usage

``` r
buoyancy_frequency(
  env_dat,
  shallow,
  deep,
  depths,
  frequency = FALSE,
  name = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- shallow:

  name of the density column for the shallower level, in kg/m^3 or as
  sigma-theta; both give the same answer since only the difference
  enters

- deep:

  name of the density column for the deeper level

- depths:

  the two depths in metres, shallower first, positive downwards

- frequency:

  return \\N\\ in s^-1 instead of \\N^2\\ in s^-2. `N` is `NA` wherever
  the column is unstable, since the root of a negative number is not a
  frequency

- name:

  name for the new column

## Value

`env_dat` with a stratification column

## Details

This is the quantity that decides how much work mixing has to do, and it
governs a great deal of what plankton do. A strongly stratified column
keeps a surface bloom in the light and cuts it off from nutrients below;
a weakly stratified one lets both mix.
[`vertical_gradient()`](https://camilleross.org/derivoce/reference/vertical_gradient.md)
approximates the same idea with a temperature difference, which is a
reasonable proxy only where salinity is uniform — and in the Gulf of
Maine it is not, because Scotian Shelf inflow is fresh enough to
stratify water that is barely warmer at the surface.

## Getting the two levels

Two routes, depending on the source.

`datamatch::accessCopernicus()` takes a `depth` argument and returns one
level per call, so a second call at a different depth gives the deeper
level and the two are joined as columns on the same points:

    surface <- datamatch::accessCopernicus(vars = c("SST", "SSS"), depth = c(0, 1), ...)
    deep    <- datamatch::accessCopernicus(vars = c("SST", "SSS"), depth = c(90, 100), ...)

`accessHYCOM()` and `accessFVCOM()` instead serve the bottom as named
variables, `BOTT` and `BOTS` beside `SST` and `SSS`, so one fetch
carries both levels and the second depth is the water depth rather than
a chosen one. That makes a surface-to-bottom **density** difference
available, where a surface-and-bottom fetch from Copernicus gives a
temperature difference only: Copernicus serves `bottomT` but no bottom
salinity.

Either way, compute a density for each level with
[`potential_density()`](https://camilleross.org/derivoce/reference/potential_density.md)
and pass both here.

## What a two-level estimate is and is not

It is a **bulk** stratification over the layer between the two depths,
not a profile. Where a sharp pycnocline sits between the levels, the
true peak \\N^2\\ at that interface is far larger than the layer average
reported here, and the depth of the peak is invisible.

The choice of levels therefore does much of the work, and the trap is
the mixed layer. If the two straddle its base, \\N^2\\ is dominated by
how much of the well-mixed layer was included rather than by the
stratification itself, and it will vary through the season for that
reason alone. Compare the levels against `MLD` before reading a seasonal
cycle into the result.

Negative values are returned rather than suppressed. They mean denser
water sits above lighter, which is genuine convective instability in
winter and otherwise usually a sign that the two levels are not what you
think.

## See also

[`potential_density()`](https://camilleross.org/derivoce/reference/potential_density.md),
[`vertical_gradient()`](https://camilleross.org/derivoce/reference/vertical_gradient.md),
[`eady_growth_rate()`](https://camilleross.org/derivoce/reference/eady_growth_rate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- potential_density(env, "SST", "SSS", name = "rho_surface")
env <- potential_density(env, "SST_deep", "SSS_deep", name = "rho_deep")
env <- buoyancy_frequency(env, "rho_surface", "rho_deep",
                          depths = c(0.494, 92.326))
} # }
```
