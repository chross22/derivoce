# Potential density of seawater

Density is what stratification, mixing and buoyancy actually depend on,
and temperature alone is a poor stand-in for it wherever salinity
varies. The Gulf of Maine is one of those places: Scotian Shelf inflow
arrives cold *and* fresh, and the two pull the density in opposite
directions, so a cold anomaly can be either denser or lighter than the
water it displaces depending on how fresh it is.

## Usage

``` r
potential_density(
  env_dat,
  temperature = "SST",
  salinity = "SSS",
  sigma = TRUE,
  name = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- temperature:

  name of the temperature column, in degrees C. Should be a potential
  temperature

- salinity:

  name of the salinity column, in PSU

- sigma:

  return sigma-theta, that is density minus 1000. `FALSE` returns
  absolute density

- name:

  name for the new column

## Value

`env_dat` with a density column, in kg/m^3

## Details

Uses the UNESCO (1983) equation of state at one atmosphere. Copernicus
`thetao` is already a potential temperature, so applying it here gives
potential density directly, conventionally reported as sigma-theta:
density in kg/m^3 minus 1000.

"One atmosphere" describes the *pressure* the density is referenced to,
not the depth the temperature and salinity came from. Potential density
is exactly the quantity you want for water sampled at depth — it is what
that water would weigh if brought to the surface, which is what makes
two levels comparable. `datamatch::accessCopernicus()` takes a `depth`
argument, so this applies to any level, and
[`buoyancy_frequency()`](https://camilleross.org/derivoce/reference/buoyancy_frequency.md)
uses two of them.

## What it is not

This is the density a parcel would have if brought to the surface. It is
the right quantity for comparing water masses and for deciding what
floats over what, and it deliberately ignores pressure, so it is **not**
in-situ density and should not be used where the compressibility of deep
water matters.

Applying it to a temperature that is not a potential temperature gives
in-situ surface density instead, which is the same number at the surface
and increasingly wrong with depth.

## Range

The polynomial is fitted over roughly -2 to 40 degrees C and 0 to 42
PSU. Outside that it still returns a number, and that number is an
extrapolation of a fit rather than a density, so values beyond the range
are warned about. Fresh water from a river mouth and ice-melt surface
layers are the usual culprits.

## References

UNESCO (1983). Algorithms for computation of fundamental properties of
seawater. *UNESCO Technical Papers in Marine Science* **44**.

## See also

[`water_mass_fraction()`](https://camilleross.org/derivoce/reference/water_mass_fraction.md),
[`vertical_gradient()`](https://camilleross.org/derivoce/reference/vertical_gradient.md),
[`cell_anomaly()`](https://camilleross.org/derivoce/reference/cell_anomaly.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- potential_density(env)

# Density anomalies say more about inflow than temperature anomalies do,
# because they combine the cold and the fresh into one number.
env <- cell_anomaly(env, "sigma_theta")
} # }
```
