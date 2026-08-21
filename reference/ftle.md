# Finite-Time Lyapunov Exponent

Measures how fast two initially neighbouring water parcels separate over
an integration window. Ridges of high FTLE mark Lagrangian coherent
structures — the transport barriers and convergence lines that organise
surface flow.

## Usage

``` r
ftle(
  env_dat,
  u = "UO",
  v = "VO",
  integration_days = 14,
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

- integration_days:

  length of the integration window, in days

- direction:

  `"backward"` (attracting structures, the default) or `"forward"`
  (repelling structures)

- step_hours:

  Runge-Kutta step size, in hours. Smaller is more accurate and slower;
  the default resolves a particle crossing roughly one grid cell per
  step at typical shelf speeds.

- name:

  name for the new column

## Value

`env_dat` with an FTLE column added, in units of 1/day. Grid-boundary
cells and particles that leave the domain or run aground are `NA`.

## Details

The direction matters, and is the reason this takes an argument rather
than picking one:

- **Backward** (the default) integrates into the past, so ridges are
  *attracting* structures: places where water arriving from different
  origins converges. These are where buoyant particles and weak swimmers
  pile up, so this is usually the direction of interest for plankton
  aggregation.

- **Forward** integrates into the future, so ridges are *repelling*
  structures: barriers that neighbouring parcels straddle and are pulled
  apart by. These separate water masses rather than concentrating them.

Particles are seeded on the covariate grid at each time step, advected
through the time-varying velocity field with a fourth-order Runge-Kutta
scheme, and the resulting flow map is differentiated to give the
Cauchy-Green strain tensor. The FTLE is then

\$\$\sigma = \frac{1}{\|T\|} \ln \sqrt{\lambda\_{max}(C)}\$\$

where \\\lambda\_{max}(C)\\ is the largest eigenvalue of the strain
tensor and \\T\\ the integration time. Units are inverse days.

## Choosing an integration time

`integration_days` sets the scale of structure resolved: too short and
no separation has accumulated, too long and the ridges blur as particles
wander far from where they started. For mesoscale ocean features, days
to a few weeks is the usual range.

Note that monthly-mean velocity fields, which is what
`datamatch::accessCopernicus()` returns for a `P1M` product, have
already averaged away the eddies that generate the sharpest structures.
FTLE computed from monthly means is a real diagnostic of the mean
circulation, but it is not the same quantity as FTLE from daily or
hourly fields, and its ridges are correspondingly smoother. Prefer a
daily (`P1D`) product where the choice exists.

## References

Haller, G. (2015). Lagrangian coherent structures. *Annual Review of
Fluid Mechanics*, 47, 137-162.

## Examples

``` r
if (FALSE) { # \dontrun{
env <- datamatch::accessCopernicus(
  dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
  vars = c("UO", "VO"), ...
)
# Where water converges - candidate aggregation sites
env <- ftle(env, integration_days = 14)       # UO and VO are the defaults
# Where water masses are pulled apart - transport barriers
env <- ftle(env, direction = "forward")

# At depth: one model level per fetch
deep <- datamatch::accessCopernicus(vars = c("UO", "VO"), depth = c(100, 100), ...)
deep <- ftle(deep, integration_days = 14)
} # }
```
