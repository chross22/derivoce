# derivoce

<!-- badges: start -->
[![R-CMD-check](https://github.com/chross22/derivoce/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chross22/derivoce/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Derived environmental covariates for ocean species distribution models — spatial
and temporal gradients, time-integrated variables, and temporal lags, computed
from gridded ocean data.

Takes the output of
[`datamatch::accessEnvDat()`](https://github.com/chross22/datamatch) (an `sf`
point object per time step) and returns the same shape with derived columns
added, so derived covariates flow into a model alongside the variables they came
from — for example into
[`taupatch`](https://github.com/chross22/taupatch).

<details>
<summary><b>Contents</b></summary>

- [Installation](#installation)
  - [datamatch, and why that line does not install it](#datamatch-and-why-that-line-does-not-install-it)
  - [Downloading data needs a Python client](#downloading-data-needs-a-python-client)
  - [Optional extras](#optional-extras)
- [Usage](#usage)
  - [Column names](#column-names)
- [What it computes](#what-it-computes)
  - [Horizontal gradients](#horizontal-gradients)
  - [Vertical gradients](#vertical-gradients)
  - [Temporal gradients, lags, and integrals](#temporal-gradients-lags-and-integrals)
  - [Fronts, contours, and flow structure](#fronts-contours-and-flow-structure)
    - [FTLE or FSLE?](#ftle-or-fsle)
- [Requirements on the input](#requirements-on-the-input)
  - [If the input was resampled](#if-the-input-was-resampled)
  - [Combinations that don't make sense are flagged](#combinations-that-dont-make-sense-are-flagged)
  - [If satellite gaps were filled](#if-satellite-gaps-were-filled)
- [Still to come](#still-to-come)

Longer form, with the reasoning behind each quantity:
[`docs/methods.md`](docs/methods.md).

</details>

## Installation

```r
# install.packages("remotes")
remotes::install_github("chross22/derivoce")
```

### datamatch, and why that line does not install it

[`datamatch`](https://github.com/chross22/datamatch) is what fetches the data
these functions derive from, but it is a **suggested** dependency rather than a
hard one, so it is not pulled in above. Install it separately — it is not on
CRAN:

```r
remotes::install_github("chross22/datamatch")
```

Or get both at once — `dependencies = TRUE` includes suggested packages, and
derivoce's `Remotes:` field tells `remotes` where datamatch lives:

```r
remotes::install_github("chross22/derivoce", dependencies = TRUE)
```

The reason it is only suggested: everything here operates on the *shape*
`accessEnvDat()` returns — an `sf` POINT object with one row per (grid point,
time step) and `YEAR`/`MONTH`/`DAY` columns — not on datamatch itself. An object
of that shape from any source works, and the test suite builds its own without
touching Copernicus. So derivoce runs without datamatch installed; you just
usually don't want it to.

One version note. The default column names here — `SST`, `BOTT`, `UO`, `VO`,
`DEPTH` — are datamatch *catalog* names, which arrived with its variable
dictionary. Against an older datamatch that returns raw Copernicus codes, the
defaults will not match and you should pass column names explicitly:

```r
eke(env, u = "uo", v = "vo")
```

If `datamatch::variable_dictionary()` exists, you are on a new enough version.

### Downloading data needs a Python client

datamatch shells out to `copernicusmarine`, the official Copernicus client. It
is not an R package and is not installed with either package:

```bash
pip install copernicusmarine
copernicusmarine login
```

`login` needs a free [Copernicus Marine
account](https://data.marine.copernicus.eu/register) and only has to be run
once. If the client is installed but R cannot find it — common when it lives in
a conda environment RStudio does not inherit the `PATH` of — point at it in
`~/.Rprofile`:

```r
options(datamatch.copernicusmarine = "~/miniconda3/bin/copernicusmarine")
```

Nothing in derivoce contacts Copernicus. This is only needed to *get* the data;
deriving from data you already have needs none of it.

### Optional extras

`distance_to_shore()` needs `rnaturalearth`, and `resolution = "large"`
additionally needs `rnaturalearthhires`, which is not on CRAN:

```r
install.packages("rnaturalearth")
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
```

Both are optional: the functions that need them say so and name the install
command if they are missing, rather than failing obscurely.

## Usage

```r
library(derivoce)

env <- datamatch::accessEnvDat(
  vars = c("SST", "BOTT"),               # product and dataset inferred
  years = 2003:2017, months = 1:12,
  bounding_box = list(xmin = -76, xmax = -65, ymin = 35, ymax = 45)
)

env <- env |>
  horizontal_gradient("SST") |>          # SST_grad, degrees C per km
  vertical_gradient() |>                 # SST - BOTT, the defaults
  temporal_gradient("SST") |>
  lag_covariate("SST") |>                # SST_lag1
  integrate_covariate("SST")             # SST_int
```

Every function takes and returns the same object, so they compose in a pipe and
nothing has to be reassembled afterwards.

### Column names

Every default here is a
[`datamatch::variable_dictionary()`](https://github.com/chross22/datamatch) name
— `SST`, `BOTT`, `UO`, `VO`, `DEPTH` — because `accessEnvDat()` now returns
columns under the names you requested rather than under Copernicus codes. So
`vertical_gradient()`, `eke()`, `current_speed()`, `ftle()`, and `fsle()` need
no column arguments on a dictionary fetch, and `distance_to_isobath()` needs
none once `datamatch::attach_bathymetry()` has added `DEPTH`:

```r
env <- datamatch::accessEnvDat(vars = c("UO", "VO"), ...) |>
  current_speed() |>                     # reads UO/VO
  eke()

# Raw Copernicus codes still work; name the columns if you used them
env <- datamatch::accessEnvDat(vars = c("uo", "vo"), ...) |>
  eke(u = "uo", v = "vo")
```

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
and near zero where the column is well mixed. Both inputs (`SST`, `BOTT`) come
from the same Copernicus dataset, so this needs no extra download.

Pass a `depth` column to get a per-metre rate instead of a total difference.
`datamatch::attach_bathymetry()` is where that column comes from:

```r
bathy <- datamatch::fetch_bathymetry(bounding_box = bb)
env   <- datamatch::attach_bathymetry(env, bathy, "DEPTH")
env   <- vertical_gradient(env, depth = "DEPTH")   # degrees C per metre
```

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
  `distance_to_isobath()` reads the `DEPTH` column that
  `datamatch::attach_bathymetry()` adds, so no separate bathymetry source is
  needed.
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
- `distance_to_shore()` — kilometres to the nearest coast, from Natural Earth.
  Static, so it is computed once per location and shared across time steps. This
  is the older pipeline's `dist`. A broad proxy for several things at once —
  depth, terrestrial input, tidal mixing, larval retention all covary with it —
  which makes it a useful covariate and a poor explanation.

#### FTLE or FSLE?

They ask inverse questions. **F**inite**T**ime**L**yapunov**E**exponents fix the integration time and measure how
far parcels separate; **F**inite**S**pace**L**yapunov**E**exponents fix a separation and measure how long it takes.

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

### If the input was resampled

`datamatch` can now put two products on one grid
(`upscale_grid()`/`downscale_grid()`) or change the time step
(`upscale_time()`/`downscale_time()`). Resampled output keeps the regular lattice
and the `YEAR`/`MONTH`/`DAY` stamping, so everything here runs on it — but two of
those directions change what a derivative means:

- **A spatial gradient of a downscaled variable measures the source grid.**
  Rendering a 0.25° field at 4 km gives it 4 km cells and no new information, so
  `horizontal_gradient()` returns the steps between the original coarse cells,
  divided by the new smaller spacing. With the default `nearest`/`step` methods
  that is visibly blocky; with `bilinear` it is a smooth field that looks like a
  measured gradient and is not. Compute gradients on the native grid and upscale
  the *result* if you need it coarser.
- **`temporal_gradient()` on time-interpolated data inherits the interpolation.**
  `downscale_time(method = "linear")` puts a constant slope between each pair of
  source steps, so the per-day rate of change is a property of the interpolant,
  not of the ocean.

Aggregating (`upscale_grid()`, `upscale_time()`) is the safe direction, and
`min_coverage` matters for `integrate_covariate()`: a partial period returned as
`NA` is excluded from the running total rather than counted as a low value.

### Combinations that don't make sense are flagged

`vars = NULL` means *every* covariate column, and `datamatch` now adds columns
that are not covariates to differentiate. Asking for a derivative that cannot
carry information gets a warning naming the column:

```r
lag_covariate(env, "DEPTH")
#> Warning: Static covariate(s) in a temporal operation: DEPTH.
#>   These hold the same value at each location in every time step, so a lag
#>   reproduces the column, a temporal gradient is zero, and an integral is a
#>   running multiple of it.
#>   Seafloor terrain from datamatch::attach_bathymetry() (DEPTH, SLOPE, ASPECT,
#>   TPI) is static in this way. Name the variables you meant, rather than
#>   letting `vars = NULL` take every column.

horizontal_gradient(env, "NAO")
#> Warning: Spatially uniform covariate(s) in a spatial operation: NAO.
#>   These take one value across the whole grid within each time step, so a
#>   horizontal gradient is zero everywhere ...
```

Two degeneracies, checked independently and only against the operation they
actually break:

| | Temporal (`lag_covariate()`, `temporal_gradient()`, `integrate_covariate()`) | Spatial (`horizontal_gradient()`, `distance_to_contour()`) |
|---|---|---|
| **Static** — `DEPTH`, `SLOPE`, `ASPECT`, `TPI` | warns | fine — this is how you get slope |
| **Spatially uniform** — `NAO`, `AO`, `AMO`, `PDO` | fine — a lagged index is a real covariate | warns |

The test is on the data, not on a list of known column names, so a variable that
happens to be constant in your particular extract is caught too, and a static
covariate from some other source is caught without the package having heard of
it.

**Non-numeric columns behave differently**, because nothing can be computed at
all rather than computed uselessly. `fill_satellite_gaps()` adds a
`<var>_source` factor recording where each value came from; naming it explicitly
is an error, while `vars = NULL` skips it silently — a caller who didn't name it
didn't mean it, and failing the whole call would make the `NULL` default
unusable on any gap-filled object.

The warnings are a safety net, not a substitute for saying what you mean. Name
your variables once the object carries anything beyond a plain `accessEnvDat()`
fetch.

### If satellite gaps were filled

Satellite and model chlorophyll differ in mean and variance, so a gradient
computed across a satellite/model boundary partly measures the change of source
rather than a feature of the ocean. `rescale = TRUE` reduces that step; the
`<var>_source` column is what tells you where the seams are.

Leaving the gaps unfilled has the opposite cost: a central difference needs both
neighbours, so every cloud hole erases a one-cell ring around itself, and
`integrate_covariate()` accumulates over whatever steps survived. Neither is
free — but the gap-filled version at least records which it is.

## Still to come

- **Vertical gradients from a depth profile** — the current one is
  surface-minus-bottom. A true `dT/dz` needs several model levels in one object,
  and `accessEnvDat()` returns one level per call: a `depth` range spanning
  several levels is now refused with a clear error rather than mislabelling
  columns. Stacking per-level fetches is the workaround, and doing it inside
  `vertical_gradient()` is the work.
- **Gulf Stream Index** — unlike NAO and AMO (now in
  [datamatch](https://github.com/chross22/datamatch), via
  `attach_climate_index()`), it has several competing definitions published in
  papers rather than at a stable URL, so it needs a decision about which one
