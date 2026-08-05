# derivoce

Derived environmental covariates for ocean species distribution models — spatial
and temporal gradients, time-integrated variables, and temporal lags, computed
from gridded ocean data.

Takes the output of
[`datamatch::accessEnvDat()`](https://github.com/chross22/datamatch) (an `sf`
point object per time step) and returns the same shape with derived columns
added, so derived covariates flow into a model alongside the variables they came
from — for example into
[`taupatch`](https://github.com/chross22/taupatch).

## Installation

```r
# install.packages("remotes")
remotes::install_github("chross22/derivoce")
```

## Usage

```r
library(derivoce)

env <- datamatch::accessEnvDat(
  product_id = "GLOBAL_MULTIYEAR_PHY_001_030",
  dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1M-m",
  vars = c("thetao", "bottomT"),
  years = 2003:2017, months = 1:12,
  bounding_box = list(xmin = -76, xmax = -65, ymin = 35, ymax = 45)
)

env <- env |>
  horizontal_gradient("thetao") |>       # thetao_grad, degrees C per km
  vertical_gradient("thetao", "bottomT") |>
  temporal_gradient("thetao") |>
  lag_covariate("thetao") |>             # thetao_lag1
  integrate_covariate("thetao")          # thetao_int
```

Every function takes and returns the same object, so they compose in a pipe and
nothing has to be reassembled afterwards.

## What it computes

### Horizontal gradients

`horizontal_gradient()` gives the magnitude of a covariate's change with
distance — what picks out fronts, the convergence zones where plankton
aggregate. Optionally the eastward and northward components too.

Results are in **covariate units per kilometre**, not per degree. That matters: a
degree of longitude is about 83 km at 42°N but 111 km at the equator, so a
per-degree gradient is stretched by latitude and not comparable across a study
area. The longitude spacing is recomputed per grid row; ignoring that introduces
roughly a 30% error at 45°.

This is a deliberate departure from the `raster::terrain()` call the older
pipeline used for `sst_grad` and `uv_grad`. `terrain()` treats its input as an
elevation in the same units as the coordinates and returns a slope *angle*, which
is dimensionally meaningless for a field measured in °C. What this returns is a
real rate of change with a real unit.

### Vertical gradients

`vertical_gradient()` is the surface-minus-bottom difference in each cell: a
stratification index, large where a warm surface layer sits over cold deep water
and near zero where the column is well mixed. Both inputs (`thetao`, `bottomT`)
come from the same Copernicus dataset, so this needs no extra download. Pass a
`depth` column to get a per-metre rate instead of a total difference.

### Temporal gradients, lags, and integrals

- `temporal_gradient()` — rate of change between consecutive time steps at each
  location, per step, per day, or per month. How fast conditions are shifting, as
  distinct from what they currently are.
- `lag_covariate()` — the value *n* steps back. Populations respond with a delay:
  a bloom feeds the animals sampled a month later, not those sampled during it.
  `n = 1` reproduces the older pipeline's `lag_sst`.
- `integrate_covariate()` — accumulation over preceding steps. A survey samples
  the food built up since the season began, not the food present that instant.
  The default `window = "year"` reproduces the older pipeline's `int_chl`
  (running sum from January, reset each year); a numeric window gives a rolling
  total that does not reset.

Locations are matched by coordinate rather than row order, so time steps need not
list their points in the same order.

### Fronts, contours, and flow structure

- `distance_to_front()` — how far each point is from the nearest front. Usually a
  better predictor than the local gradient: a station in a smooth patch has zero
  gradient whether the nearest front is 2 km or 200 km away, and those are very
  different places to be.
- `distance_to_contour()` / `distance_to_isobath()` — distance to where a
  covariate crosses a value. Plankton track the shelf break, and "20 km inshore of
  the 100 m isobath" locates that better than "depth = 85 m".
- `ftle()` / `fsle()` — Lyapunov exponents, forward or backward. **Backward** (the
  default) finds *attracting* structures, where water converges and material
  accumulates; **forward** finds *repelling* structures, the transport barriers.
  Depth-resolved versions work today by fetching velocities at a chosen depth —
  Copernicus GLORYS carries `uo`/`vo` on 50 levels.
- `eke()` — eddy kinetic energy, with an explicit choice of what the anomaly is
  measured against (record mean, monthly climatology, or a rolling window),
  because that choice decides what counts as an eddy rather than as mean flow.
- `current_speed()` — speed from u/v, reproducing the older pipeline's `uv`; pipe
  it through `horizontal_gradient()` for `uv_grad`.

#### FTLE or FSLE?

They ask inverse questions. **FTLE** fixes the integration time and measures how
far parcels separate; **FSLE** fixes a separation and measures how long it takes.

```r
env <- ftle(env, integration_days = 14)   # "how much separation in 14 days?"
env <- fsle(env, final_separation = 50)   # "how long to separate by 50 km?"
```

Prefer **FSLE** when the *spatial* scale is what matters, or when the domain spans
very different flow speeds. An FTLE map with one fixed integration time resolves
fine structure where the flow is fast and only coarse structure where it is slow,
so across a shelf with an energetic break current and a sluggish interior, ridge
intensity partly encodes current speed rather than frontal activity — and a model
cannot separate the two. FSLE asks the same question everywhere.

Prefer **FTLE** when the *timescale* is what matters and can be named — a
retention time, a cohort's accumulation window, time since a bloom. FSLE has
nowhere to encode that.

Two caveats outweigh the choice: monthly fields have already averaged away the
eddies that make sharp structures, and plankton are not passive surface tracers.
[`docs/methods.md`](docs/methods.md) covers both, along with how each quantity is
computed, which choices were deliberate, and where results differ from the older
pipeline.

## Requirements on the input

Spatial derivatives are only defined on a grid, so `horizontal_gradient()`
requires points on a **regular lon/lat lattice** — which gridded products
(Copernicus and similar) are, and scattered observations are not. Irregular input
is rejected rather than interpolated, since silently gridding it would produce a
gradient field that looks plausible and is mostly interpolation artifact.

Central differences are undefined at the grid edge, so boundary cells come back
`NA`. So do the first *n* time steps of a lag, and the first step of a temporal
gradient.

## Still to come

- **FSLE** — the finite-*size* counterpart to FTLE, resolving structures at a
  chosen spatial scale rather than a chosen timescale. The advection machinery is
  reusable; only the stopping criterion differs.
- **Distance to shore** — static, needs a coastline source
- **Climate indices** (NAO, AMO, Gulf Stream Index) — a retrieval problem rather
  than a derivation, so probably better placed in `datamatch`
