# Eddy kinetic energy

The kinetic energy carried by the *departure* of the flow from its mean:

## Usage

``` r
eke(env_dat, u = "UO", v = "VO", reference = "record", name = "EKE")
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- u:

  name of the eastward velocity column, in m/s

- v:

  name of the northward velocity column, in m/s

- reference:

  `"record"`, `"climatology"`, or a positive integer number of time
  steps for a centred rolling mean

- name:

  name for the new column

## Value

`env_dat` with an EKE column added, in m^2/s^2

## Details

\$\$EKE = \tfrac{1}{2}(u'^2 + v'^2), \quad u' = u - \bar{u}\$\$

High EKE marks an energetic, variable region — meanders, rings, and
eddies — as distinct from a strong but steady current, which carries
high total kinetic energy but little eddy energy.

## What the anomaly is measured against

EKE is only defined relative to a mean, and the choice of mean decides
what counts as an "eddy". This is a scientific decision, not an
implementation detail, so it is an argument with no silently-correct
default:

- `"record"` (the default) subtracts the mean over the whole time series
  at each location. This is the textbook definition. Because the
  seasonal cycle of the mean circulation is left in the anomaly, a
  region whose currents merely strengthen every summer will register as
  energetic.

- `"climatology"` subtracts a separate mean for each calendar month, so
  the repeatable seasonal cycle is removed and only departures *from the
  usual conditions for that month* count. Use this when the question is
  about anomalous years rather than about which places are energetic.

- An integer subtracts a centred rolling mean over that many time steps,
  so slow drift in the mean state is removed and only variability faster
  than the window survives.

Means are computed per location, so a spatially varying mean circulation
is removed correctly rather than being smeared into the anomaly.

## Examples

``` r
if (FALSE) { # \dontrun{
env <- eke(env)                             # against the record mean
env <- eke(env, reference = "climatology")  # seasonal cycle removed
env <- eke(env, reference = 12)             # 12-step rolling mean

# If the fetch used Copernicus codes rather than catalog names
env <- eke(env, u = "uo", v = "vo")
} # }
```
