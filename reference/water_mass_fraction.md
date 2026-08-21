# Water-mass fraction from two endmembers

Estimates what fraction of each cell is the first of two named water
masses, by projecting its temperature and salinity onto the mixing line
between the two endmembers in T-S space.

## Usage

``` r
water_mass_fraction(
  env_dat,
  endmembers,
  temperature = "SST",
  salinity = "SSS",
  residual = FALSE,
  name = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- endmembers:

  a list of exactly two named `c(temperature, salinity)` vectors. The
  fraction returned is of the **first**.

- temperature:

  name of the temperature column, in degrees C

- salinity:

  name of the salinity column, in PSU

- residual:

  also add the distance from the mixing line, as `<name>_residual`

- name:

  name for the fraction column. `NULL` uses the first endmember's name,
  lowercased, with `_frac`.

## Value

`env_dat` with a fraction column in `[0, 1]`, and optionally a residual
column

## Details

The Gulf of Maine is the standard setting for this. Its deep water is a
mixture of Warm Slope Water and Labrador Slope Water in proportions that
vary episodically, and the proportion sets the nutrient supply: Labrador
Slope Water is colder, fresher, and nutrient-poor, Warm Slope Water
warmer, saltier, and nutrient-rich (Townsend et al. 2015). Surface water
is instead a mixture of Scotian Shelf Water with slope water, which is
the same calculation with different endmembers.

## How the fraction is computed

A cell's `(T, S)` is projected onto the line joining the two endmembers,
and the fraction is where it lands: 1 at the first endmember, 0 at the
second.

Temperature and salinity are on different scales, so each axis is
normalised by the endmember separation in that variable before
projecting. Without that the axis with the larger numerical spread would
dominate the answer for reasons of units rather than of oceanography.

Values are clamped to `[0, 1]`. A cell beyond an endmember is either a
third water mass or an endmember choice that does not bracket the data,
and reporting a fraction of 1.4 would dress up that problem as a
measurement.

## Check the residual

Projection always returns a fraction, even for water that is not a
mixture of these two masses at all. `residual = TRUE` adds the distance
from the mixing line in the same normalised units, which is the number
that says whether the fraction means anything: near zero, the cell sits
on the line and the mixture is a good description; large, and it does
not.

This is why the endmembers are an argument with no default. They vary by
region, season, and year, and a wrong pair produces confident nonsense.

## References

Townsend DW, Pettigrew NR, Thomas MA, Neary MG, McGillicuddy DJ,
O'Donnell J (2015). Water masses and nutrient sources to the Gulf of
Maine. *Journal of Marine Research* **73**, 93-122.

## See also

[`derived_indices()`](https://camilleross.org/derivoce/reference/derived_indices.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Deep water: which slope water is filling the basins?
env <- water_mass_fraction(env, endmembers = list(
  LSW = c(temperature = 6,  salinity = 34.4),
  WSW = c(temperature = 12, salinity = 35.4)
), residual = TRUE)

# Surface water: how much of this is Scotian Shelf Water?
env <- water_mass_fraction(env, endmembers = list(
  SSW   = c(temperature = 2,  salinity = 32.0),
  slope = c(temperature = 10, salinity = 35.0)
))
} # }
```
