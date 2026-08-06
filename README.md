# derivoce

<!-- badges: start -->
[![R-CMD-check](https://github.com/chross22/derivoce/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chross22/derivoce/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Derived oceanographic covariates for species distribution models: spatial and
temporal gradients, time-integrated variables, temporal lags, and fluid dynamics,
computed from gridded ocean data.

It takes the output of
[`datamatch::accessEnvDat()`](https://github.com/chross22/datamatch), an `sf`
point object per time step. It returns the same shape, with derived columns
added. 

<details>
<summary><b>Contents</b></summary>

- [Installation](#installation)
  - [datamatch, and why that line does not install it](#datamatch-and-why-that-line-does-not-install-it)
  - [Downloading data needs a Python client](#downloading-data-needs-a-python-client)
  - [Optional extras](#optional-extras)
- [Usage](#usage)
- [What it computes](#what-it-computes)
  - [Horizontal gradients](#horizontal-gradients)
  - [Vertical gradients](#vertical-gradients)
  - [Temporal gradients, lags, and integrals](#temporal-gradients-lags-and-integrals)
  - [Fronts, contours, and flow structure](#fronts-contours-and-flow-structure)
    - [FTLE or FSLE?](#ftle-or-fsle)
- [Column names](#column-names)
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

[`datamatch`](https://github.com/chross22/datamatch) fetches the data these
functions derive from. It is a **suggested** dependency rather than a hard one,
so the line above does not install it. It is not on CRAN either, so get it from
GitHub:

```r
remotes::install_github("chross22/datamatch")
```

Or install both at once. `dependencies = TRUE` includes suggested packages, and
derivoce's `Remotes:` field tells `remotes` where to find datamatch:

```r
remotes::install_github("chross22/derivoce", dependencies = TRUE)
```

Why is it only suggested? Everything here works on the *shape* that
`accessEnvDat()` returns. That shape is an `sf` POINT object with one row per
grid point and time step, plus `YEAR`, `MONTH`, and `DAY` columns. Nothing here
calls datamatch itself. An object of that shape works no matter where it came
from, and the test suite builds its own without touching Copernicus. So derivoce
runs fine without datamatch installed. You just usually don't want it to.

### Downloading data needs a Python client

datamatch downloads through `copernicusmarine`, the official Copernicus client.
It is not an R package, and needs to be installed via command line, like below:

```bash
pip install copernicusmarine
copernicusmarine login
```

`login` needs a free [Copernicus Marine
account](https://data.marine.copernicus.eu/register), and you only have to run
it once. Sometimes the client installs fine but R cannot find it. That is common
when it lives in a conda environment whose `PATH` RStudio does not inherit. Point
at it directly in `~/.Rprofile`:

```r
options(datamatch.copernicusmarine = "~/miniconda3/bin/copernicusmarine")
```

Nothing in derivoce contacts Copernicus. You only need this to *get* the data.
Deriving from data you already have needs none of it.

### Optional extras

`distance_to_shore()` needs `rnaturalearth`. The `resolution = "large"` option
also needs `rnaturalearthhires`, which is not on CRAN:

```r
install.packages("rnaturalearth")
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
```

Both are optional. If one is missing, the function that needs it says so and
gives you the install command.

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

Every function takes and returns the same object. So they compose in a pipe, and
nothing has to be reassembled afterwards.

Notice that `vertical_gradient()` above needed no arguments. [Column
names](#column-names) explains why, once the functions it applies to have been
introduced.

## What it computes

### Horizontal gradients

`horizontal_gradient()` gives the magnitude of a covariate's change with
distance. That is what picks out fronts, the convergence zones where plankton
aggregate. It can return the eastward and northward components too.

Results are in **covariate units per kilometre**, not per degree. That matters.
A degree of longitude is about 83 km at 42°N, but 111 km at the equator. A
per-degree gradient is therefore stretched by latitude, and not comparable across
a study area. The longitude spacing is recomputed for every grid row. Ignoring
that introduces roughly a 30% error at 45°.

This is a deliberate departure from the `raster::terrain()` call the older
pipeline used for `sst_grad` and `uv_grad`. `terrain()` treats its input as an
elevation measured in the same units as the coordinates, and returns a slope
*angle*. That is dimensionally meaningless for a field measured in °C. What this
returns is a real rate of change, with a real unit.

### Vertical gradients

`vertical_gradient()` is the surface-minus-bottom difference in each cell, which
works as a stratification index. It is large where a warm surface layer sits over
cold deep water, and near zero where the column is well mixed. Both inputs
(`SST` and `BOTT`) come from the same Copernicus dataset, so this needs no extra
download.

Pass a `depth` column to get a per-metre rate instead of a total difference.
`datamatch::attach_bathymetry()` is where that column comes from:

```r
bathy <- datamatch::fetch_bathymetry(bounding_box = bb)
env   <- datamatch::attach_bathymetry(env, bathy, "DEPTH")
env   <- vertical_gradient(env, depth = "DEPTH")   # degrees C per metre
```

### Temporal gradients, lags, and integrals

- `temporal_gradient()` gives the rate of change between consecutive time steps
  at each location, per step, per day, or per month. It measures how fast
  conditions are shifting, as distinct from what they currently are.
- `lag_covariate()` gives the value *n* steps back. Populations respond with a
  delay. A bloom feeds the animals sampled a month later, not those sampled
  during it. `n = 1` reproduces the older pipeline's `lag_sst`.
- `integrate_covariate()` accumulates over preceding steps. A survey samples the
  food built up since the season began, not the food present at that instant.
  The default `window = "year"` reproduces the older pipeline's `int_chl`, a
  running sum from January that resets each year. A numeric window gives a
  rolling total that never resets.

Locations are matched by coordinate rather than by row order. So time steps need
not list their points in the same order.

### Fronts, contours, and flow structure

- `distance_to_front()` measures how far each point is from the nearest front.
  It is usually a better predictor than the local gradient. A station in a smooth
  patch has zero gradient whether the nearest front is 2 km away or 200 km away,
  and those are very different places to be.
- `distance_to_contour()` and `distance_to_isobath()` measure the distance to
  where a covariate crosses a value. Plankton track the shelf break, and "20 km
  inshore of the 100 m isobath" locates that better than "depth = 85 m".
  `distance_to_isobath()` reads the `DEPTH` column that
  `datamatch::attach_bathymetry()` adds, so you need no separate bathymetry
  source.
- `ftle()` and `fsle()` compute Lyapunov exponents, forward or backward.
  **Backward** is the default. It finds *attracting* structures, where water
  converges and material accumulates. **Forward** finds *repelling* structures,
  the transport barriers. Depth-resolved versions work today: fetch velocities at
  a chosen depth, since Copernicus GLORYS carries `uo` and `vo` on 50 levels.
- `eke()` computes eddy kinetic energy. You choose explicitly what the anomaly is
  measured against, whether the record mean, a monthly climatology, or a rolling
  window. That choice decides what counts as an eddy rather than as mean flow.
- `current_speed()` gives speed from u and v, reproducing the older pipeline's
  `uv`. Pipe it through `horizontal_gradient()` for `uv_grad`.
- `distance_to_shore()` gives kilometres to the nearest coast, from Natural
  Earth. It is static, so it is computed once per location and shared across time
  steps. This is the older pipeline's `dist`. It is a broad proxy for several
  things at once. Depth, terrestrial input, tidal mixing, and larval retention
  all covary with it. That makes it a useful covariate and a poor explanation.

#### Why FTLE and FSLE come back mostly NA

This is the most common surprise, and it is not a failure. Both work by
following parcels through the velocity field. A parcel that reaches the edge of
the data has no velocity to follow, so it is retired and its cell returns `NA`.

That costs a margin around the domain of roughly **speed × integration time**. A
typical shelf speed of 0.15 m/s over the default 14 days is about 180 km. On a
500 km box that removes a third of the field. On a 1° box it removes all of it:

```r
ftle(env, integration_days = 14)
#> Warning: ftle() returned no values at all: every point is NA.
#>   A 14-day integration at this field's median speed (0.2 m/s) carries a
#>   parcel about 250 km, and the domain is 82 by 110 km. 507 of 507 particles
#>   left the velocity field before the window was up.
#>   Shorten integration_days, or fetch a larger bounding box. Either way a
#>   margin of roughly speed x integration_days is lost along the upstream
#>   edge, so the usable area is always smaller than the area fetched.
```

So **fetch a bounding box larger than your study area**, by about that margin.
Backward integration loses the upstream edge, forward the downstream one.

FSLE has a second way to return nothing: parcels that stay in the domain but
never separate by `final_separation`. A near-uniform flow carries pairs along
together rather than pulling them apart, and monthly means have already averaged
away the eddies that do the pulling. Its warning distinguishes the two causes,
because they call for opposite fixes — a shorter `max_days` for parcels lost to
the edge, a longer one for parcels that never separated.

The warning fires only when almost everything is `NA`. Losing a margin is normal
and is not worth interrupting for.

#### FTLE or FSLE?

They ask inverse questions. **F**inite-**T**ime **L**yapunov **E**xponents fix
the integration time and measure how far parcels separate. **F**inite-**S**ize
**L**yapunov **E**xponents fix a separation and measure how long it takes.

```r
env <- ftle(env, integration_days = 14)   # "how much separation in 14 days?"
env <- fsle(env, final_separation = 50)   # "how long to separate by 50 km?"
```

Prefer **FSLE** when the *spatial* scale is what matters, or when the domain
spans very different flow speeds. An FTLE map uses one fixed integration time. It
resolves fine structure where the flow is fast, and only coarse structure where
it is slow. Across a shelf with an energetic break current and a sluggish
interior, ridge intensity then partly encodes current speed rather than frontal
activity. A model cannot separate the two. FSLE asks the same question
everywhere.

Prefer **FTLE** when the *timescale* is what matters and can be named: a
retention time, a cohort's accumulation window, or time since a bloom. FSLE has
nowhere to encode that.

Two caveats outweigh the choice. Monthly fields have already averaged away the
eddies that make sharp structures. And plankton are not passive surface tracers.
[`docs/methods.md`](docs/methods.md) covers both. It also covers how each
quantity is computed, which choices were deliberate, and where results differ
from the older pipeline.

## Column names

Every default above is a datamatch catalog name: `SST` and `BOTT` for
`vertical_gradient()`, `UO` and `VO` for `eke()`, `current_speed()`, `ftle()`,
and `fsle()`, and `DEPTH` for `distance_to_isobath()`.

That is not a coincidence. `accessEnvDat()` returns columns under the names you
asked for, rather than under Copernicus codes. So a dictionary fetch needs no
column arguments at all.

Here is the same `eke()` call twice. Only the fetch differs:

```r
bb <- list(xmin = -76, xmax = -65, ymin = 35, ymax = 45)

# Asked for catalog names, so the columns are UO and VO.
# That is exactly what eke() looks for by default.
env <- datamatch::accessEnvDat(vars = c("UO", "VO"), years = 2010, months = 1:12,
                               bounding_box = bb)
env <- eke(env)

# Asked for Copernicus codes, so the columns are uo and vo. Same data and
# same function, but now the columns have to be named.
env <- datamatch::accessEnvDat(vars = c("uo", "vo"), years = 2010, months = 1:12,
                               bounding_box = bb)
env <- eke(env, u = "uo", v = "vo")
```

`distance_to_isobath()` is the one exception. Its `DEPTH` column does not come
from `accessEnvDat()` at all, but from `datamatch::attach_bathymetry()`.

## Requirements on the input

Spatial derivatives are only defined on a grid. So `horizontal_gradient()`
requires points on a **regular lon/lat lattice**. Gridded products like
Copernicus are regular, and scattered observations are not. Irregular input is
rejected rather than interpolated. Silently gridding it would produce a gradient
field that looks plausible but is mostly interpolation artifact.

Central differences are undefined at the grid edge, so boundary cells come back
`NA`. So do the first *n* time steps of a lag, and the first step of a temporal
gradient.

### If the input was resampled

`datamatch` can now put two products on one grid, with `upscale_grid()` and
`downscale_grid()`. It can also change the time step, with `upscale_time()` and
`downscale_time()`. Resampled output keeps the regular lattice and the
`YEAR`/`MONTH`/`DAY` stamping, so everything here runs on it. But two of those
directions change what a derivative means.

- **A spatial gradient of a downscaled variable measures the source grid.**
  Rendering a 0.25° field at 4 km gives it 4 km cells and no new information.
  `horizontal_gradient()` then returns the steps between the original coarse
  cells, divided by the new smaller spacing. With the default `nearest` and
  `step` methods that is visibly blocky. With `bilinear` it is a smooth field
  that looks like a measured gradient but is not. Compute gradients on the native
  grid, then upscale the *result* if you need it coarser.
- **`temporal_gradient()` on time-interpolated data inherits the
  interpolation.** `downscale_time(method = "linear")` puts a constant slope
  between each pair of source steps. The per-day rate of change is then a
  property of the interpolant, not of the ocean.

Aggregating with `upscale_grid()` or `upscale_time()` is the safe direction. One
interaction is worth knowing. `min_coverage` matters for
`integrate_covariate()`: a partial period returned as `NA` is excluded from the
running total, rather than counted as a low value.

### Combinations that don't make sense are flagged

`vars = NULL` means *every* covariate column. And `datamatch` now adds columns
that are not covariates to differentiate. If you ask for a derivative that cannot
carry information, you get a warning naming the column:

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

There are two degeneracies. Each is checked independently, and only against the
operation it actually breaks:

| | Temporal (`lag_covariate()`, `temporal_gradient()`, `integrate_covariate()`) | Spatial (`horizontal_gradient()`, `distance_to_contour()`) |
|---|---|---|
| **Static**: `DEPTH`, `SLOPE`, `ASPECT`, `TPI` | warns | fine, this is how you get slope |
| **Spatially uniform**: `NAO`, `AO`, `AMO`, `PDO`, `LCR` | fine, a lagged index is real | warns |

The test looks at the data, not at a list of known column names. So a variable
that happens to be constant in your particular extract is caught too. And a
static covariate from some other source is caught without the package ever having
heard of it.

**Non-numeric columns behave differently.** With these, nothing can be computed
at all, rather than computed uselessly. `fill_satellite_gaps()` adds a
`<var>_source` factor recording where each value came from. Naming it explicitly
is an error. But `vars = NULL` skips it silently. If you did not name it, you did
not mean it, and failing the whole call would make the `NULL` default unusable on
any gap-filled object.

The warnings are a safety net, not a substitute for saying what you mean. Name
your variables once the object carries anything beyond a plain `accessEnvDat()`
fetch.

### If satellite gaps were filled

Satellite and model chlorophyll differ in mean and variance. So a gradient
computed across a satellite/model boundary partly measures the change of source,
rather than a feature of the ocean. `rescale = TRUE` reduces that step, and the
`<var>_source` column tells you where the seams are.

Leaving the gaps unfilled has the opposite cost. A central difference needs both
neighbours, so every cloud hole erases a one-cell ring around itself. And
`integrate_covariate()` accumulates over whatever steps survived. Neither option
is free. But the gap-filled version at least records which one you got.

## Still to come

- **Vertical gradients from a depth profile.** The current one is
  surface-minus-bottom. A true `dT/dz` needs several model levels in one object.
  `accessEnvDat()` returns one level per call, and a `depth` range spanning
  several levels is now refused with a clear error rather than mislabelling
  columns. Stacking per-level fetches is the workaround. Doing that inside
  `vertical_gradient()` is the work.
- **Gulf Stream Index.** NAO and AMO are now in
  [datamatch](https://github.com/chross22/datamatch), via
  `attach_climate_index()`. The Gulf Stream Index is harder. It has several
  competing definitions, published in papers rather than at a stable URL. So it
  needs a decision about which one to use.
