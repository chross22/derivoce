# Finite-Size Lyapunov Exponent

Measures how *quickly* neighbouring water parcels reach a chosen
separation, rather than how far apart they get in a chosen time. Where
[`ftle()`](https://camilleross.org/derivoce/reference/ftle.md) fixes the
clock and measures distance, FSLE fixes the distance and measures the
clock:

## Usage

``` r
fsle(
  env_dat,
  u = "UO",
  v = "VO",
  final_separation = 50,
  initial_separation = NULL,
  max_days = 60,
  direction = c("backward", "forward"),
  step_hours = 6,
  name = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return, on a regular lon/lat grid,
  containing eastward and northward velocity columns

- u:

  name of the eastward velocity column, in m/s

- v:

  name of the northward velocity column, in m/s

- final_separation:

  target separation \\\delta_f\\, in km. This is the scale-selectivity
  knob: it sets the size of structure resolved.

- initial_separation:

  starting separation \\\delta_0\\, in km; defaults to roughly one grid
  cell

- max_days:

  give up after this long. Must exceed the time a typical parcel pair
  needs, or most of the field returns `NA`.

- direction:

  `"backward"` (attracting structures, the default) or `"forward"`
  (repelling structures)

- step_hours:

  Runge-Kutta step size, in hours

- name:

  name for the new column

## Value

`env_dat` with an FSLE column added, in 1/day. Parcels that never reach
`final_separation` within `max_days` are `NA`.

## Details

\$\$\lambda = \frac{\ln(\delta_f / \delta_0)}{\tau}\$\$

with \\\delta_0\\ the starting separation, \\\delta_f\\ the target, and
\\\tau\\ the time taken to reach it. Units are 1/day.

## Why choose this over FTLE

The difference is **scale selectivity**, and it matters most in a domain
whose flow speed varies a lot from place to place.

An FTLE map with a fixed integration time resolves fine structure where
the flow is fast and only coarse structure where it is slow. Across a
shelf with an energetic break current and a sluggish interior, ridge
intensity then partly encodes background current speed rather than
frontal activity, and a model cannot tell the two apart.

FSLE asks the same question everywhere — "how long to separate by
\\\delta_f\\?" — so results are comparable between energetic and quiet
regions, and \\\delta_f\\ can be set to a scale that means something
biologically: a patch size, a predator's search radius, a survey's
resolution.

The cost is that FSLE has no natural place to encode a biologically
meaningful *timescale*. When the relevant duration is known — a
retention time, a cohort's accumulation window —
[`ftle()`](https://camilleross.org/derivoce/reference/ftle.md) with that
integration time is the better tool.

## How it is computed

Each grid point is seeded with two companion particles, offset east and
north by `initial_separation`. All three are advected together, and the
separation of each pair is checked after every step. \\\tau\\ is the
first time either pair reaches `final_separation`.

Parcels that never separate that far within `max_days` return `NA`: the
question "how long to reach this separation" simply has no answer for
them, and substituting `max_days` would report a slow separation rate
where there was none at all.

## References

d'Ovidio, F., Fernandez, V., Hernandez-Garcia, E., & Lopez, C. (2004).
Mixing structures in the Mediterranean Sea from finite-size Lyapunov
exponents. *Geophysical Research Letters*, 31(17).

## See also

[`ftle()`](https://camilleross.org/derivoce/reference/ftle.md), and
`docs/methods.md` for when to prefer which

## Examples

``` r
if (FALSE) { # \dontrun{
# Structures at the 50 km scale, comparable across the whole shelf
env <- datamatch::accessCopernicus(vars = c("UO", "VO"), ...)
env <- fsle(env, final_separation = 50)       # UO and VO are the defaults

# A scale matched to a predator's search radius
env <- fsle(env, final_separation = 10, max_days = 60)
} # }
```
