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
  - [Installing datamatch](#installing-datamatch)
  - [The Copernicus client](#the-copernicus-client)
  - [Optional extras](#optional-extras)
- [Usage](#usage)
- [What it computes](#what-it-computes)
  - [Horizontal gradients](#horizontal-gradients)
  - [Vertical gradients](#vertical-gradients)
  - [Temporal gradients, lags, and integrals](#temporal-gradients-lags-and-integrals)
    - [Lag by calendar time, not by position](#lag-by-calendar-time-not-by-position)
  - [Fronts, contours, and flow structure](#fronts-contours-and-flow-structure)
  - [FTLE or FSLE?](#ftle-or-fsle)
  - [Regional indices](#regional-indices)
    - [Three ways to measure the same inflow](#three-ways-to-measure-the-same-inflow)
    - [Using your own line or box](#using-your-own-line-or-box)
    - [Before you use these](#before-you-use-these)
    - [A note on "Follows"](#a-note-on-follows)
- [Column names](#column-names)
- [Requirements on the input](#requirements-on-the-input)
- [Warnings you may see](#warnings-you-may-see)
  - [Mostly NA from FTLE or FSLE](#mostly-na-from-ftle-or-fsle)
  - [A derivative that cannot carry information](#a-derivative-that-cannot-carry-information)
- [Resampled and gap-filled input](#resampled-and-gap-filled-input)
- [Still to come](#still-to-come)
- [References](#references)
  - [Data sources](#data-sources)
  - [Software](#software)
  - [Keeping these current](#keeping-these-current)
  - [Citing derivoce](#citing-derivoce)

Longer form, with the reasoning behind each quantity:
[`docs/methods.md`](docs/methods.md).

</details>

## Installation

```r
# install.packages("remotes")
remotes::install_github("chross22/derivoce")
```

### Installing datamatch

[`datamatch`](https://github.com/chross22/datamatch) fetches the data these
functions derive from. It is a **suggested** dependency, so the line above does
not install it. It is not on CRAN either:

```r
remotes::install_github("chross22/datamatch")

# or get both at once
remotes::install_github("chross22/derivoce", dependencies = TRUE)
```

Why only suggested? Everything here works on the *shape* `accessEnvDat()`
returns: an `sf` POINT object with one row per grid point and time step, plus
`YEAR`, `MONTH`, and `DAY`. Nothing calls datamatch itself. So an object of that
shape works whatever produced it, and the tests build their own.

### The Copernicus client

datamatch downloads through `copernicusmarine`, the official Copernicus client.
It is not an R package, and needs to be installed from the command line:

```bash
pip install copernicusmarine
copernicusmarine login
```

`login` needs a free [Copernicus Marine
account](https://data.marine.copernicus.eu/register) and is a one-off.

Sometimes R cannot find the client afterwards. That is common when it sits in a
conda environment whose `PATH` RStudio does not inherit. Point at it directly in
`~/.Rprofile`:

```r
options(datamatch.copernicusmarine = "~/miniconda3/bin/copernicusmarine")
```

Nothing in derivoce contacts Copernicus. This is only for *getting* the data.

### Optional extras

`distance_to_shore()` needs `rnaturalearth`, and `resolution = "large"` also
needs `rnaturalearthhires`, which is not on CRAN:

```r
install.packages("rnaturalearth")
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
```

If either is missing, the function that needs it says so and gives the install
command.

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

Every function takes and returns the same object, so they compose in a pipe.

`vertical_gradient()` needed no arguments there. [Column
names](#column-names) explains why.

## What it computes

### Horizontal gradients

`horizontal_gradient()` gives the magnitude of a covariate's change with
distance. That is what picks out fronts, the convergence zones where plankton
aggregate. It can return the eastward and northward components too.

Results are in **covariate units per kilometre**, not per degree. That matters:
a degree of longitude is about 83 km at 42°N but 111 km at the equator, so a
per-degree gradient is stretched by latitude and not comparable across a study
area. The longitude spacing is recomputed for every grid row.

This is a deliberate departure from the `raster::terrain()` call the original
pipeline used for `sst_grad` and `uv_grad` (Ross et al. 2023). `terrain()` returns a slope *angle*,
which is dimensionally meaningless for a field in °C. This returns a real rate of
change with a real unit.

### Vertical gradients

`vertical_gradient()` is the difference between two temperature columns in each
cell, a stratification index. It is large where a warm surface layer sits over
cold deep water, and near zero where the column is mixed. It defaults to
surface-minus-bottom, and both of those come from the same Copernicus dataset,
so the usual case needs no extra download.

The two columns are arguments rather than fixed, so any two levels work:
`accessEnvDat()` takes a `depth` and returns one level per call, so a second
fetch gives a temperature at whatever depth you choose. Where you know the two
depths, `buoyancy_frequency()` is the better measure — a temperature difference
stands in for stratification only where salinity is uniform, and Scotian Shelf
inflow is fresh enough to stratify water barely warmer at the surface.

Pass a `depth` column for a per-metre rate instead of a total difference.
`datamatch::attach_bathymetry()` supplies that column:

```r
bathy <- datamatch::fetch_bathymetry(bounding_box = bb)
env   <- datamatch::attach_bathymetry(env, bathy, "DEPTH")
env   <- vertical_gradient(env, depth = "DEPTH")   # degrees C per metre
```

### Temporal gradients, lags, and integrals

- `temporal_gradient()` gives the rate of change between consecutive steps, per
  step, per day, or per month. How fast conditions are shifting, as distinct from
  what they are.
- `lag_covariate()` gives the value *n* steps back. Populations respond with a
  delay: a bloom feeds the animals sampled a month later, not those sampled
  during it. Ross et al. (2023) used a one-month SST lag, which is
  `by = "month"` here.
- `integrate_covariate()` accumulates over preceding steps. A survey samples the
  food built up since the season began, not the food present at that instant. The
  default `window = "year"` reproduces the `int_chl` of Ross et al. (2023),
  chlorophyll integrated from January and reset each year. A numeric window rolls
  without resetting.

Locations are matched by coordinate, not row order, so time steps need not list
their points in the same order.

#### Lag by calendar time, not by position

`lag_covariate()` counts steps by default, which is only unambiguous when the
series is evenly spaced and complete. `by` counts calendar time instead:

```r
lag_covariate(env, "CHL", n = 3, by = "month")   # CHL_lag3month
lag_covariate(env, "SST", n = 1, by = "year")    # same month, last year
lag_covariate(env, "SST", n = 30, by = "day")    # daily products
```

Prefer a calendar unit whenever the lag means something biological. "Three
months ago" is a claim about the organism. "Three steps ago" is a claim about
how the data was fetched, and the two stop agreeing the moment a month is
missing. In a monthly series missing April, `by = "step"` makes March the
predecessor of May, so a one-step lag is quietly a two-month one. `by = "month"`
returns `NA` there instead.

`n` may be a vector, which is what an autoregressive design needs:

```r
lag_covariate(env, "SST", n = 1:3, by = "year")
# adds SST_lag1year, SST_lag2year, SST_lag3year
```

With `by = "year"` this holds the calendar month fixed and varies only the year,
so the seasonal cycle drops out and what remains is interannual.

### Anomalies, extremes, and density

- `cell_anomaly()` removes each cell's own mean, so what is left is the
  departure rather than the geography. Eight degrees is cold for the southern
  Gulf of Maine and warm for the Scotian Shelf, and a model given raw
  temperature has to learn that before it can use the difference.
  `standardize = TRUE` divides by the cell's own variability too, which makes
  departures comparable across the domain at the cost of the magnitude.
  `detrend = TRUE` also removes the long-term trend, which matters in a warming
  shelf sea: an anomaly that still contains the trend largely encodes *which
  year it is*. The cell-wise counterpart of `box_anomaly()`.
- `decompose_covariate()` splits a series into its parts — a trend, a repeating
  seasonal cycle, and the residual — additively, so each can be used or
  inspected on its own. `slope` gives the rate of change per year, one number
  per cell. The trend and the cycle are fitted **together**, because removed one
  after the other each absorbs part of the other.
- `index_series()` collapses a broadcast per-step column back to one row per
  time step. The region-scale indices below all compute one value per step and
  repeat it on every row so the object keeps its shape; this pulls the series
  back out for plotting or export. It refuses to collapse a column that varies
  across the grid, because that would keep one arbitrary cell and discard the
  map.
- `marine_heatwave()` flags periods unusually warm — or, with
  `direction = "cold"`, unusually cold — for the time of year, after Hobday et
  al. (2016, 2018), with intensity, duration, cumulative intensity and the four
  categories. Whether an animal was there in an ordinary summer or a heatwave
  summer is often more useful than the temperature itself.
- `potential_density()` gives sigma-theta from temperature and salinity, by the
  UNESCO (1983) equation of state. Density, not temperature, is what
  stratification and buoyancy depend on, and the two part company exactly where
  it matters here: Scotian Shelf inflow is cold *and* fresh, which pull density
  in opposite directions.
- `buoyancy_frequency()` gives N² between two depths — the real stratification
  measure, where `vertical_gradient()` only approximates it with a temperature
  difference. Needs a second `accessEnvDat()` call at a deeper `depth`, joined
  as columns.
- `eady_growth_rate()` says how fast baroclinic instability grows, from the
  vertical shear, the stratification and the Coriolis parameter. High where a
  sheared, weakly stratified flow can overturn — the shelf-break front and the
  edges of warm-core rings — so it complements `detect_eddies()`: that finds
  eddies which exist, this finds where conditions favour making them. (Eady is
  a person, not a spelling of "eddy".)

### Fronts, contours, and flow structure

- `distance_to_front()` measures how far each point is from the nearest front,
  usually a better predictor than the local gradient. Fronts are found by
  thresholding the gradient, after Belkin and O'Reilly (2009). A station in a smooth patch
  has zero gradient whether the nearest front is 2 km or 200 km away, and those
  are very different places to be.
- `distance_to_contour()` and `distance_to_isobath()` measure distance to where a
  covariate crosses a value. Plankton track the shelf break, and "20 km inshore of
  the 100 m isobath" locates that better than "depth = 85 m".
- `ftle()` and `fsle()` compute Lyapunov exponents (Haller 2015; d'Ovidio et al.
  2004). **Backward** (the default)
  finds *attracting* structures, where water converges and material accumulates.
  **Forward** finds *repelling* structures, the transport barriers. For
  depth-resolved versions, fetch velocities at a chosen depth: Copernicus GLORYS
  carries `uo` and `vo` on 50 levels.
- `front_frequency()` asks how *reliably* a front sits in a place, rather than
  how far one was at a moment. Fronts move: a cell frontal in one step of twenty
  caught a passing filament, while one frontal in fifteen sits on a shelf-break
  or tidal mixing front, and only the second aggregates plankton reliably enough
  for a predator to learn it.
- `flow_deformation()` gives vorticity, divergence, the strain components,
  the Okubo-Weiss parameter and the Rossby number from the velocity gradients.
  Instantaneous and local, so it sits between `eke()`, which needs a series, and
  the Lyapunov exponents, which need trajectories: separating an eddy interior
  from the filaments around it costs one pass rather than an integration.
- `detect_eddies()` goes from a field to objects: it groups the connected cells
  where rotation beats strain into individual eddies and describes each one, so
  a cell carries not "how eddy-like is the flow here" but "you are inside an
  eddy, it turns this way, and it is this big". Polarity is the part that earns
  its keep — cyclonic cores upwell and often concentrate plankton, anticyclonic
  ones downwell, and a covariate that only says "eddy" averages the two together
  and can easily find nothing.
- `distance_to_eddy()` completes the `distance_to_*` family. A station 5 km
  outside a rotating core and one 300 km away are different places, and an
  inside/outside flag scores both zero. `polarity` narrows it to cyclonic or
  anticyclonic.
- `residence_time()` releases a particle at every point in a box and measures
  how long it stays. Long residence means a retentive place where anything with
  a life stage measured in weeks can complete it. Read the censoring note in
  `?residence_time` before averaging the result.
- `eke()` computes eddy kinetic energy. You choose what the anomaly is measured
  against: the record mean, a monthly climatology, or a rolling window. That
  choice decides what counts as an eddy rather than mean flow.
- `current_speed()` gives speed from u and v, the `uv` of Ross et al. (2023).
  Their `uv_grad` is the spatial derivative of that speed field, so it is
  `current_speed()` then `horizontal_gradient()` on the result. Differentiating
  `u` and `v` separately and combining afterwards is a different quantity.
- `distance_to_shore()` gives kilometres to the nearest coast, from [Natural
  Earth](https://www.naturalearthdata.com/). Static, so it is computed once per location and shared across time steps.
  A broad proxy for several things at once: depth, terrestrial input, tidal
  mixing, and larval retention all covary with it. Useful as a covariate, poor as
  an explanation.

### FTLE or FSLE?

They ask inverse questions. **F**inite-**T**ime **L**yapunov **E**xponents fix
the integration time and measure how far parcels separate. **F**inite-**S**ize
**L**yapunov **E**xponents fix a separation and measure how long it takes.

```r
env <- ftle(env, integration_days = 14)   # "how much separation in 14 days?"
env <- fsle(env, final_separation = 50)   # "how long to separate by 50 km?"
```

Prefer **FSLE** when the *spatial* scale is what matters, or when the domain
spans very different flow speeds. One fixed integration time resolves fine
structure where the flow is fast and coarse structure where it is slow, so ridge
intensity ends up partly encoding current speed rather than frontal activity.
FSLE asks the same question everywhere.

Prefer **FTLE** when the *timescale* is what matters and can be named: a
retention time, a cohort's accumulation window, time since a bloom. FSLE has
nowhere to put that.

Two caveats outweigh the choice. Monthly fields have already averaged away the
eddies that make sharp structures, and plankton are not passive surface tracers.
[`docs/methods.md`](docs/methods.md) covers both.

### Regional indices

Most functions here give you a value for every grid cell. These four give you one
number per month for a whole region, like a climate index. They answer "how much
water came in this month", not "what was it like here".

Two currents feed the Gulf of Maine, and they carry very different water:

![Two doorways into the Gulf of Maine](docs/figures/water-mass-origins.png)

Cold, fresh, nutrient-poor water rounds Cape Sable from the Scotian Shelf. Warm,
salty, nutrient-rich water comes in deep through the Northeast Channel. The two
take turns, so which one is dominant changes what the Gulf is like that season.
That is why they are two indices and not one.

```r
env <- scotian_shelf_inflow(env)       # m^2/s, positive = into the Gulf
env <- northeast_channel_inflow(env)
```

#### Three ways to measure the same inflow

You can ask three different questions about Scotian Shelf water arriving, and the
literature asks all three. They are not interchangeable:

| Question | Function | Needs | Follows |
|---|---|---|---|
| How much water crossed this line? | `scotian_shelf_inflow()` | `UO`, `VO` | Feng et al. 2016; Wang et al. 2022 |
| How much of the water here came from there? | `water_mass_fraction()` | `SST`, `SSS` | Townsend et al. 2015 |
| Did the water here get fresher? | `eastern_gom_salinity()` | `SSS` | Grodsky et al. 2025 |

The first measures the flow itself, and is the only one that gives you a
direction. The second measures what is present rather than what moved, which is
what matters for nutrients, and it works on data with no currents in it. The
third is the simplest and the least specific: it tells you conditions changed,
not that water moved. Freshening could equally be rain or runoff.

Use more than one and disagreement is informative. Strong inflow with no
freshening means the water that arrived was not unusually fresh, which tells you
something about the Scotian Shelf that year.

```r
env <- water_mass_fraction(env, endmembers = list(
  LSW = c(temperature = 6,  salinity = 34.4),
  WSW = c(temperature = 12, salinity = 35.4)
), residual = TRUE)

env <- eastern_gom_salinity(env)
```

`derived_indices()` lists all of them with their sources.
`derived_indices(markdown = TRUE)` gives you a table to paste elsewhere.

#### Using your own line or box

The named indices have fixed geometry, because an index named after a place is
defined by that place. For anywhere else, use the general versions:

```r
env <- section_transport(env, from = c(-66.5, 43.3), to = c(-65.6, 42.6))
env <- box_anomaly(env, "SSS", box = list(xmin = -68, xmax = -66,
                                          ymin = 43, ymax = 44.5))
```

#### Before you use these

**They flip sign in summer.** Positive through winter, negative from June to
September, at both sections. That is the real surface circulation, not a bug.
Treat them as winter indices.

**The numbers are not comparable to published transports.** These integrate one
surface layer along a line. A mooring array integrates the full depth of the
section, so the figures differ by orders of magnitude. Read these as "more or
less than usual", not as a flux.

**The Northeast Channel changed after 2000.** Gulf Stream warm-core rings drive
slope water in, and ring formation nearly doubled around then (Silver et al.
2023). A record spanning 2000 covers two different regimes, so check any
long-term relationship on each side separately.

**Check the residual on `water_mass_fraction()`.** It always returns a fraction,
even for water that is not a mix of your two endmembers at all. `residual = TRUE`
is how you find out whether the answer means anything.

We chose the section endpoints ourselves by testing them against real currents;
they are not from any paper.
[`docs/methods.md`](docs/methods.md#how-the-named-sections-were-placed) shows how,
and [`docs/section-placement-diagnostics.R`](docs/section-placement-diagnostics.R)
re-runs the test on your own data.

#### A note on "Follows"

It means we implemented the idea, not that we reproduce the published series.
Each function computes from whatever data you give it, so the numbers will differ
from the paper's. Cite the paper for the concept and describe your own inputs.
Sources are also available as `as.data.frame(derived_indices())$source`, and all
work cited anywhere here is listed under [References](#references) at the end.

## Column names

Every default above is a datamatch catalog name: `SST` and `BOTT` for
`vertical_gradient()`, `UO` and `VO` for `eke()`, `current_speed()`, `ftle()`,
and `fsle()`, and `DEPTH` for `distance_to_isobath()`.

`accessEnvDat()` returns columns under the names you asked for, rather than under
Copernicus codes, so a dictionary fetch needs no column arguments. Here is the
same `eke()` call twice, differing only in the fetch:

```r
bb <- list(xmin = -76, xmax = -65, ymin = 35, ymax = 45)

# Catalog names, so the columns are UO and VO - what eke() expects by default.
env <- datamatch::accessEnvDat(vars = c("UO", "VO"), years = 2010, months = 1:12,
                               bounding_box = bb)
env <- eke(env)

# Copernicus codes, so the columns are uo and vo and have to be named.
env <- datamatch::accessEnvDat(vars = c("uo", "vo"), years = 2010, months = 1:12,
                               bounding_box = bb)
env <- eke(env, u = "uo", v = "vo")
```

`distance_to_isobath()` is the exception: its `DEPTH` column comes from
`datamatch::attach_bathymetry()`, not from `accessEnvDat()`.

## Requirements on the input

Spatial derivatives are only defined on a grid, so `horizontal_gradient()`
requires points on a **regular lon/lat lattice**. Gridded products are regular.
Scattered observations are not. Irregular input is rejected rather than interpolated,
because a gradient computed from interpolated data mostly measures the
interpolation.

Central differences are undefined at the grid edge, so boundary cells come back
`NA`. So do the first *n* steps of a lag and the first step of a temporal
gradient. These are genuine absences, not failures.

## Warnings you may see

### Mostly NA from FTLE or FSLE

This is the most common surprise, and it is not a failure. Both follow parcels
through the velocity field, and a parcel that reaches the edge of the data has no
velocity left to follow, so its cell returns `NA`.

That costs a margin of roughly **speed × integration time** around the domain. At
a shelf speed of 0.15 m/s the default 14 days is about 180 km, which removes a
third of a 500 km box and all of a 1° one:

```r
ftle(env, integration_days = 14)
#> Warning: ftle() returned no values at all: every point is NA.
#>   A 14-day integration at this field's median speed (0.2 m/s) carries a
#>   parcel about 250 km, and the domain is 82 by 110 km. 507 of 507 particles
#>   left the velocity field before the window was up.
#>   Shorten integration_days, or fetch a larger bounding box...
```

So **fetch a bounding box larger than your study area**, by about that margin.
Backward integration loses the upstream edge, forward the downstream one.

FSLE can also return nothing for a second reason: parcels that stay in the domain
but never separate by `final_separation`. Its warning tells the two apart,
because they need opposite fixes. Parcels lost to the edge want a shorter
`max_days`. Parcels that never separated want a longer one.

The warning only fires when almost everything is `NA`. Losing a margin is normal.

### A derivative that cannot carry information

`vars = NULL` means *every* covariate column, and datamatch now attaches columns
that are not covariates to differentiate. Asking for a derivative that cannot say
anything gets a warning naming the column:

```r
lag_covariate(env, "DEPTH")
#> Warning: Static covariate(s) in a temporal operation: DEPTH.
#>   These hold the same value at each location in every time step, so a lag
#>   reproduces the column, a temporal gradient is zero, and an integral is a
#>   running multiple of it...
```

Two degeneracies, each checked only against the operation it actually breaks:

| | Temporal (`lag_covariate()`, `temporal_gradient()`, `integrate_covariate()`) | Spatial (`horizontal_gradient()`, `distance_to_contour()`) |
|---|---|---|
| **Static**: `DEPTH`, `SLOPE`, `ASPECT`, `TPI` | warns | fine, this is how you get slope |
| **Spatially uniform**: `NAO`, `AO`, `AMO`, `PDO`, `LCR`, `AMOC` | fine, a lagged index is real | warns |

The test looks at the data, not at a list of known names, so a variable that
happens to be constant in your extract is caught too.

**Non-numeric columns are an error instead**, since nothing can be computed at
all. `fill_satellite_gaps()` adds a `<var>_source` factor. Naming it explicitly
fails, while `vars = NULL` skips it silently.

These warnings are a safety net, not a substitute for naming your variables once
the object carries more than a plain `accessEnvDat()` fetch.

## Resampled and gap-filled input

datamatch can put two products on one grid (`upscale_grid()`,
`downscale_grid()`) or change the time step (`upscale_time()`,
`downscale_time()`). Resampled output keeps the regular lattice and the
`YEAR`/`MONTH`/`DAY` stamping, so everything here runs on it. Two directions
change what a derivative means:

- **A spatial gradient of a downscaled variable measures the source grid.**
  Rendering a 0.25° field at 4 km adds cells, not information, so the gradient is
  the step between the original coarse cells divided by the new smaller spacing.
  Derive on the native grid and upscale the *result* instead.
- **`temporal_gradient()` on time-interpolated data measures the interpolant.**
  `downscale_time(method = "linear")` puts a constant slope between source steps,
  and that slope is what you get back.

Aggregating is the safe direction. Note that `min_coverage` interacts with
`integrate_covariate()`: a partial period returned as `NA` drops out of the
running total rather than counting as a low value.

For gap-filled satellite data, satellite and model chlorophyll differ in mean and
variance, so a gradient across a seam partly measures the change of source.
`rescale = TRUE` reduces the step, and `<var>_source` says where the seams are.
Leaving gaps unfilled costs the other way: a central difference needs both
neighbours, so every cloud hole erases a ring around itself.

## Still to come

- **Vertical gradients from a full depth profile.** The two-level case is
  covered — `vertical_gradient()` takes any two temperature columns and
  `buoyancy_frequency()` gives N² between any two depths. A true profile,
  and the potential energy anomaly that needs one, is still assembly work:
  `accessEnvDat()` returns one level per call, so a profile is several fetches
  joined as columns. Doing that inside the functions is what remains.
- **Extending the LCR index past 2014.** Tried twice and shelved. Monthly
  Copernicus fields fail for a physical reason: the averaging removes the narrow
  Labrador Current jet and the Grand Banks bifurcation the index depends on.
  Daily fields with OceanParcels clear that obstacle and still do not reproduce
  the published series — and it is not the arrival regions, the particle count
  or the domain, all of which were tested and eliminated. What is left is that
  this counts particles entering a region while the paper counts them crossing a
  hydrographic section, and those coordinates are not published.
  [`docs/lcr-extension-experiment.md`](docs/lcr-extension-experiment.md) records
  both attempts.
- **Gulf Stream Index.** NAO, AO, AMO, PDO, LCR, and AMOC are all in
  [datamatch](https://github.com/chross22/datamatch) via
  `attach_climate_index()`. The Gulf Stream Index is harder: it has several
  competing definitions published in papers rather than at a stable URL, so it
  needs a decision about which one.

## References

Work cited anywhere above, alphabetical. Each function's own `?help` carries the
references relevant to it, and `as.data.frame(derived_indices())$source` gives
them for the regional indices at runtime.

Where a function "follows" a paper, it implements that paper's idea and computes
it from whatever data you supply. None reproduces a published time series, so
cite the paper for the concept and describe your own inputs.

- Belkin IM, O'Reilly JE (2009). An algorithm for oceanic front detection in
  chlorophyll and SST satellite imagery. *Journal of Marine Systems* **78**(3),
  319–326.
  [doi:10.1016/j.jmarsys.2008.11.018](https://doi.org/10.1016/j.jmarsys.2008.11.018)
- d'Ovidio F, Fernández V, Hernández-García E, López C (2004). Mixing structures
  in the Mediterranean Sea from finite-size Lyapunov exponents. *Geophysical
  Research Letters* **31**(17).
  [doi:10.1029/2004GL020328](https://doi.org/10.1029/2004GL020328)
- Du J, Zhang WG, Li Y (2022). Impact of Gulf Stream warm-core rings on slope
  water intrusion into the Gulf of Maine. *Journal of Physical Oceanography*
  **52**(8).
  [doi:10.1175/JPO-D-21-0288.1](https://doi.org/10.1175/JPO-D-21-0288.1)
- Eady ET (1949). Long waves and cyclone waves. *Tellus* **1**(3), 33–52.
  [doi:10.3402/tellusa.v1i3.8507](https://doi.org/10.3402/tellusa.v1i3.8507)
- Feng H, Vandemark D, Wilkin J (2016). Gulf of Maine salinity variation and its
  correlation with upstream Scotian Shelf currents at seasonal and interannual
  time scales. *Journal of Geophysical Research: Oceans* **121**.
  [doi:10.1002/2016JC012337](https://doi.org/10.1002/2016JC012337)
- Grodsky SA, Vandemark D, Levin J (2025). An eastern Gulf of Maine salinity
  index for monitoring winter Scotian Shelf inflow and its relation to coastal
  and interior pathways. *Journal of Geophysical Research: Oceans* **130**(5).
  [doi:10.1029/2024JC021891](https://doi.org/10.1029/2024JC021891)
- Haller G (2015). Lagrangian coherent structures. *Annual Review of Fluid
  Mechanics* **47**, 137–162.
  [doi:10.1146/annurev-fluid-010313-141322](https://doi.org/10.1146/annurev-fluid-010313-141322)
- Hobday AJ, Alexander LV, Perkins SE, Smale DA, Straub SC, Oliver ECJ,
  Benthuysen JA, Burrows MT, Donat MG, Feng M, Holbrook NJ, Moore PJ, Scannell
  HA, Sen Gupta A, Wernberg T (2016). A hierarchical approach to defining marine
  heatwaves. *Progress in Oceanography* **141**, 227–238.
  [doi:10.1016/j.pocean.2015.12.014](https://doi.org/10.1016/j.pocean.2015.12.014)
- Hobday AJ, Oliver ECJ, Sen Gupta A, Benthuysen JA, Burrows MT, Donat MG,
  Holbrook NJ, Moore PJ, Thomsen MS, Wernberg T, Smale DA (2018). Categorizing
  and naming marine heatwaves. *Oceanography* **31**(2), 162–173.
  [doi:10.5670/oceanog.2018.205](https://doi.org/10.5670/oceanog.2018.205)
- Isern-Fontanet J, García-Ladona E, Font J (2003). Identification of marine
  eddies from altimetric maps. *Journal of Atmospheric and Oceanic Technology*
  **20**(5), 772–778.
  [doi:10.1175/1520-0426(2003)20<772:IOMEFA>2.0.CO;2](https://doi.org/10.1175/1520-0426(2003)20%3C772:IOMEFA%3E2.0.CO;2)
- Lindzen RS, Farrell B (1980). A simple approximate result for the maximum
  growth rate of baroclinic instabilities. *Journal of the Atmospheric Sciences*
  **37**(7), 1648–1654.
  [doi:10.1175/1520-0469(1980)037<1648:ASARFT>2.0.CO;2](https://doi.org/10.1175/1520-0469%281980%29037%3C1648:ASARFT%3E2.0.CO;2)
- Okubo A (1970). Horizontal dispersion of floatable particles in the vicinity of
  velocity singularities such as convergences. *Deep-Sea Research and
  Oceanographic Abstracts* **17**(3), 445–454.
  [doi:10.1016/0011-7471(70)90059-8](https://doi.org/10.1016/0011-7471(70)90059-8)
- Ramp SR, Schlitz RJ, Wright WR (1985). The deep flow through the Northeast
  Channel, Gulf of Maine. *Journal of Physical Oceanography* **15**(12),
  1790–1808.
- Ross C, Runge J, Roberts J, Brady D, Tupper B, Record N (2023). Estimating
  North Atlantic right whale prey based on *Calanus finmarchicus* thresholds.
  *Marine Ecology Progress Series* **703**, 1–16.
  [doi:10.3354/meps14204](https://doi.org/10.3354/meps14204)
- Silver A, Gangopadhyay A, Gawarkiewicz G, Fratantoni P, Clark J (2023).
  Increased Gulf Stream warm core ring formations contributes to an observed
  increase in salinity maximum intrusions on the Northeast Shelf. *Scientific
  Reports* **13**, 7538.
  [doi:10.1038/s41598-023-34494-0](https://doi.org/10.1038/s41598-023-34494-0)
- Townsend DW, Pettigrew NR, Thomas MA, Neary MG, McGillicuddy DJ, O'Donnell J
  (2015). Water masses and nutrient sources to the Gulf of Maine. *Journal of
  Marine Research* **73**, 93–122.
- UNESCO (1983). Algorithms for computation of fundamental properties of
  seawater. *UNESCO Technical Papers in Marine Science* **44**.
  [unesdoc.unesco.org](https://unesdoc.unesco.org/ark:/48223/pf0000059832)
- Wang et al. (2022). Freshwater transport in the Scotian Shelf and its impacts
  on the Gulf of Maine salinity. *Journal of Geophysical Research: Oceans*
  **127**.
  [doi:10.1029/2021JC017663](https://doi.org/10.1029/2021JC017663)
- Weiss J (1991). The dynamics of enstrophy transfer in two-dimensional
  hydrodynamics. *Physica D: Nonlinear Phenomena* **48**(2–3), 273–294.
  [doi:10.1016/0167-2789(91)90088-Q](https://doi.org/10.1016/0167-2789(91)90088-Q)

### Data sources

- **Copernicus Marine Service** supplies the gridded fields these covariates are
  derived from, chiefly the GLORYS12V1 global ocean reanalysis
  (`GLOBAL_MULTIYEAR_PHY_001_030`). Copernicus asks that products be credited in
  any publication using them; see
  <https://marine.copernicus.eu/> for the current wording and the DOI of the
  specific product and version you fetched. `datamatch::index_dictionary()` and
  `datamatch::variable_dictionary()` report which product each variable came
  from.
- **Natural Earth** provides the coastlines behind `distance_to_shore()`. Public
  domain, via `rnaturalearth`. <https://www.naturalearthdata.com/>
- **NOAA ETOPO**, via `marmap`, is the source of the depth grid used to place and
  check the named sections, and of `DEPTH` when it comes from
  `datamatch::fetch_bathymetry()`.

### Software

These do the geometric and raster work, and are worth citing alongside this
package. `citation("sf")` and so on give the current form.

- **sf** — Pebesma E (2018). Simple features for R: standardized support for
  spatial vector data. *The R Journal* **10**(1), 439–446.
  [doi:10.32614/RJ-2018-009](https://doi.org/10.32614/RJ-2018-009)
- **terra** — Hijmans R. *terra: Spatial Data Analysis*. R package.
  <https://CRAN.R-project.org/package=terra>
- **rnaturalearth** — Massicotte P, South A. *rnaturalearth: World Map Data from
  Natural Earth*. R package.
  <https://CRAN.R-project.org/package=rnaturalearth>

Version and year are deliberately omitted for the two R packages: both move with
every release, so `citation("terra")` is the answer rather than anything written
down here.

### Keeping these current

A scheduled workflow re-checks the citations each quarter: that every DOI is
still registered, that everything cited in the code or docs appears in the list
below, and that nothing in the list is cited nowhere. It opens an issue when
something needs a look, and does not try to fix anything itself, since choosing
the right replacement reference is a judgement rather than a lookup.

Run it yourself with:

```bash
Rscript inst/scripts/check_citations.R
```

It checks whether doi.org has the DOI registered, and deliberately stops there
rather than following through to the publisher. Publishers routinely answer a
scripted request with 403, and treating that as a dead reference would file a
false alarm every quarter.

### Citing derivoce

`citation("derivoce")` gives the current form. If a specific covariate follows a
published method, cite that paper too: the list above says which, and each
function's `?help` repeats it.
