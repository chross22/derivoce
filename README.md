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
    - [References](#references)
- [Column names](#column-names)
- [Requirements on the input](#requirements-on-the-input)
- [Warnings you may see](#warnings-you-may-see)
  - [Mostly NA from FTLE or FSLE](#mostly-na-from-ftle-or-fsle)
  - [A derivative that cannot carry information](#a-derivative-that-cannot-carry-information)
- [Resampled and gap-filled input](#resampled-and-gap-filled-input)
- [Still to come](#still-to-come)

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

This is a deliberate departure from the `raster::terrain()` call the older
pipeline used for `sst_grad` and `uv_grad`. `terrain()` returns a slope *angle*,
which is dimensionally meaningless for a field in °C. This returns a real rate of
change with a real unit.

### Vertical gradients

`vertical_gradient()` is the surface-minus-bottom difference in each cell, a
stratification index. It is large where a warm surface layer sits over cold deep
water, and near zero where the column is mixed. Both inputs come from the same
Copernicus dataset, so it needs no extra download.

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
  during it. `n = 1` reproduces the older pipeline's `lag_sst`.
- `integrate_covariate()` accumulates over preceding steps. A survey samples the
  food built up since the season began, not the food present at that instant. The
  default `window = "year"` reproduces `int_chl`, a running sum from January that
  resets each year. A numeric window rolls without resetting.

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

### Fronts, contours, and flow structure

- `distance_to_front()` measures how far each point is from the nearest front,
  usually a better predictor than the local gradient. A station in a smooth patch
  has zero gradient whether the nearest front is 2 km or 200 km away, and those
  are very different places to be.
- `distance_to_contour()` and `distance_to_isobath()` measure distance to where a
  covariate crosses a value. Plankton track the shelf break, and "20 km inshore of
  the 100 m isobath" locates that better than "depth = 85 m".
- `ftle()` and `fsle()` compute Lyapunov exponents. **Backward** (the default)
  finds *attracting* structures, where water converges and material accumulates.
  **Forward** finds *repelling* structures, the transport barriers. For
  depth-resolved versions, fetch velocities at a chosen depth: Copernicus GLORYS
  carries `uo` and `vo` on 50 levels.
- `eke()` computes eddy kinetic energy. You choose what the anomaly is measured
  against: the record mean, a monthly climatology, or a rolling window. That
  choice decides what counts as an eddy rather than mean flow.
- `current_speed()` gives speed from u and v, the older pipeline's `uv`. Pipe it
  through `horizontal_gradient()` for `uv_grad`.
- `distance_to_shore()` gives kilometres to the nearest coast, from Natural
  Earth. Static, so it is computed once per location and shared across time steps.
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

Some quantities describe a region rather than a cell: transport across a
section, how much of the water came from somewhere, whether a basin was fresher
than usual. These return **one value per time step**, broadcast to every row, so
they behave like a climate index rather than a map.

```r
derived_indices()                      # the catalogue, with sources
cat(derived_indices(markdown = TRUE))  # as a markdown table
```

Scotian Shelf inflow to the Gulf of Maine appears three times, because the
literature measures it three ways and they answer different questions:

| Index | Method | Needs | Measures | Follows |
|---|---|---|---|---|
| `scotian_shelf_inflow()` | transport across a fixed line off Cape Sable | `UO`, `VO` | the crossing itself, with a direction | Feng et al. 2016; Wang et al. 2022 |
| `water_mass_fraction()` | T-S endmember mixing | `SST`, `SSS` | how much of the water present came from there | Townsend et al. 2015 |
| `eastern_gom_salinity()` | box salinity anomaly | `SSS` | the most visible consequence, freshening | Grodsky et al. 2025 |
| `northeast_channel_inflow()` | transport across the Northeast Channel | `UO`, `VO` | slope water entering by the deep route | Ramp et al. 1985; Du et al. 2022; Silver et al. 2023 |

A transport is the only one that gives a flux. A water-mass fraction is what
matters for nutrients and works where velocities do not. A box anomaly is the
most robust and the least specific: it says conditions changed, not that water
moved. They disagree in useful ways, since a strong inflow with a normal
salinity anomaly means the arriving water was not unusually fresh.

`northeast_channel_inflow()` is a separate index, not a variant. The Channel is
the deep route by which slope water enters the Gulf, and it alternates
episodically with the shallow, fresh Cape Sable inflow. The contrast is the
point.

**"Follows" means the approach, not the series.** Each function follows the idea in
the cited work but computes it from whatever gridded field you supply, so none
reproduces a published time series and none should be compared to one value for
value. `eastern_gom_salinity()` is the clearest case: Grodsky et al. built their
index from SMAP satellite salinity, and this will give you something different
from a model reanalysis. Cite the paper for the concept, and describe your own
inputs.

The section endpoints are ours, not theirs. No published section is being
reproduced, and how they were arrived at is in
[`docs/methods.md`](docs/methods.md#how-the-named-sections-were-placed).

The named sections are fixed rather than arguments, since an index named for a
place is defined by that place. `section_transport()` is the general function,
and `scotian_shelf_inflow_section()` reports the geometry so it can be plotted
or checked.

Those endpoints were measured against 60 months of real GLORYS velocities rather
than chosen from a map. Around 80% of the local flow crosses each section, against
65% and 27% for an earlier pair picked by eye, the second of which ran nearly
along the channel instead of across it.

Both transports **reverse sign in summer**: positive into the Gulf through
winter, negative from June to September, at both sections independently. That is
the surface circulation, not a bug, so read these as seasonal indices and not as
a year-round inflow.
[`docs/methods.md`](docs/methods.md#how-the-named-sections-were-placed) records
how, with figures, and
[`docs/section-placement-diagnostics.R`](docs/section-placement-diagnostics.R)
re-runs the check on your own extract — worth doing for a different season or
region.

```r
env <- scotian_shelf_inflow(env)         # m^2/s, positive into the Gulf
env <- northeast_channel_inflow(env)
env <- water_mass_fraction(env, endmembers = list(
  LSW = c(temperature = 6,  salinity = 34.4),
  WSW = c(temperature = 12, salinity = 35.4)
), residual = TRUE)

# The general forms, for any line or box
env <- section_transport(env, from = c(-66.5, 43.3), to = c(-65.6, 42.6))
env <- box_anomaly(env, "SSS", box = eastern_gom_box())
```

Two cautions. A surface velocity field integrated along a line is a **proxy
for** depth-integrated transport, not a measurement of it, and the Northeast
Channel in particular is baroclinic enough that the surface can run opposite to
the deep flow. And `water_mass_fraction()` always returns a fraction, even for
water that is not a mixture of those two masses at all, so check the residual.

**The Northeast Channel inflow regime is not stationary.** Slope water entering
there is modulated by Gulf Stream warm-core rings (Du et al. 2022), and the
forcing itself changed: ring formation nearly doubled after 2000, from about 18
a year to 33, and salinity-maximum intrusions onto the Northeast Shelf
quadrupled, 72% of them coinciding with a ring offshore (Silver et al. 2023).

A long record of this index therefore spans two regimes. For anything
interannual, check whether a relationship holds before and after 2000
separately rather than assuming it is stable. It also changes what a high value
means: the same circulation can carry more slope water simply because more rings
are present.

#### References

Also available at runtime, as
`as.data.frame(derived_indices())$source`, and in each function's `?help`.

- Feng H, Vandemark D, Wilkin J (2016). Gulf of Maine salinity variation and its
  correlation with upstream Scotian Shelf currents at seasonal and interannual
  time scales. *Journal of Geophysical Research: Oceans* **121**.
  [doi:10.1002/2016JC012337](https://doi.org/10.1002/2016JC012337)
- Du J, Zhang WG, Li Y (2022). Impact of Gulf Stream warm-core rings on slope
  water intrusion into the Gulf of Maine. *Journal of Physical Oceanography*
  **52**(8). [doi:10.1175/JPO-D-21-0288.1](https://doi.org/10.1175/JPO-D-21-0288.1)
- Grodsky SA, Vandemark D, Levin J (2025). An eastern Gulf of Maine salinity
  index for monitoring winter Scotian Shelf inflow and its relation to coastal
  and interior pathways. *Journal of Geophysical Research: Oceans* **130**(5).
  [doi:10.1029/2024JC021891](https://doi.org/10.1029/2024JC021891)
- Ramp SR, Schlitz RJ, Wright WR (1985). The deep flow through the Northeast
  Channel, Gulf of Maine. *Journal of Physical Oceanography* **15**(12),
  1790–1808.
- Silver A, Gangopadhyay A, Gawarkiewicz G, Fratantoni P, Clark J (2023).
  Increased Gulf Stream warm core ring formations contributes to an observed
  increase in salinity maximum intrusions on the Northeast Shelf. *Scientific
  Reports* **13**, 7538.
  [doi:10.1038/s41598-023-34494-0](https://doi.org/10.1038/s41598-023-34494-0)
- Townsend DW, Pettigrew NR, Thomas MA, Neary MG, McGillicuddy DJ, O'Donnell J
  (2015). Water masses and nutrient sources to the Gulf of Maine. *Journal of
  Marine Research* **73**, 93–122.
- Wang et al. (2022). Freshwater transport in the Scotian Shelf and its impacts
  on the Gulf of Maine salinity. *Journal of Geophysical Research: Oceans*
  **127**. [doi:10.1029/2021JC017663](https://doi.org/10.1029/2021JC017663)

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

- **Vertical gradients from a depth profile.** The current one is
  surface-minus-bottom. A true `dT/dz` needs several model levels in one object,
  and `accessEnvDat()` returns one level per call. Stacking per-level fetches is
  the workaround. Doing that inside `vertical_gradient()` is the work.
- **Gulf Stream Index.** NAO, AO, AMO, PDO, LCR, and AMOC are all in
  [datamatch](https://github.com/chross22/datamatch) via
  `attach_climate_index()`. The Gulf Stream Index is harder: it has several
  competing definitions published in papers rather than at a stable URL, so it
  needs a decision about which one.
