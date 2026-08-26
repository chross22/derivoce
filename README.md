# derivoce

<!-- badges: start -->
[![R-CMD-check](https://github.com/chross22/derivoce/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chross22/derivoce/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Derived oceanographic covariates: spatial and temporal gradients,
time-integrated variables, temporal lags, and fluid dynamics, computed from
gridded ocean data.

It takes an `sf` point object with one row per location and time step — the
shape every [`datamatch`](https://github.com/chross22/datamatch) access function
returns, whichever source it read. It returns the same shape, with derived
columns added, so every function composes in a pipe.

Nothing here is tied to one product. datamatch currently serves seven sources —
Copernicus Marine, HYCOM, CCMP winds, FVCOM, ERDDAP, seafloor terrain and
climate indices — and derivoce works on the shape rather than the source, so an
object of that shape works whatever produced it. The one distinction that
matters is grid geometry, not provenance: see [Requirements on the
input](#requirements-on-the-input).

## Installation

```r
# install.packages("remotes")
remotes::install_github("chross22/derivoce")
```

**[`datamatch`](https://github.com/chross22/datamatch) is not a declared
dependency**, so nothing installs it for you, and it is not on CRAN:

```r
remotes::install_github("chross22/datamatch")
```

Why not even suggested? Everything here works on the *shape* datamatch's access
functions return: an `sf` POINT object with one row per location and time step,
plus `YEAR`, `MONTH`, and `DAY`. That shape is the same from every source. Nothing in this package calls datamatch — it appears
only in documentation, in `\dontrun{}` examples, and inside warning messages that
name which datamatch function a column came from. An object of that shape works
whatever produced it, and the tests build their own. Declaring it anyway had a
cost: a suggested package that resolves from no repository makes `R CMD check`
report a dependency it cannot find, for a package that is never loaded.

**The Copernicus client.** datamatch downloads through `copernicusmarine`, which
is not an R package. Nothing in derivoce contacts Copernicus — this is only for
*getting* the data:

```bash
pip install copernicusmarine
copernicusmarine login
```

If R cannot find it afterwards — common when it sits in a conda environment whose
`PATH` RStudio does not inherit — point at it directly in `~/.Rprofile`:

```r
options(datamatch.copernicusmarine = "~/miniconda3/bin/copernicusmarine")
```

**Optional extras.** `distance_to_shore()` needs `rnaturalearth`, and
`resolution = "large"` also needs `rnaturalearthhires`, which is not on CRAN. If
either is missing, the function that needs it says so and gives the install
command.

```r
install.packages("rnaturalearth")
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
```

## Usage

```r
library(derivoce)

env <- datamatch::accessCopernicus(
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

`vertical_gradient()` needed no arguments there — see [Column
names](#column-names) for why.

## What it computes

| | |
|---|---|
| **Gradients** | `horizontal_gradient()` (per **kilometre**, not per degree), `vertical_gradient()` (stratification between two levels), `temporal_gradient()` |
| **Time** | `lag_covariate()`, `integrate_covariate()`, `rolling_covariate()`, `decompose_covariate()` (trend, cycle, residual, fitted together), `index_series()` |
| **Anomalies and extremes** | `cell_anomaly()`, `box_anomaly()`, `marine_heatwave()` (Hobday categories, warm or cold) |
| **Density and stratification** | `potential_density()` (sigma-theta, UNESCO 1983), `buoyancy_frequency()` (N², the real measure), `eady_growth_rate()`, `water_mass_fraction()` |
| **Fronts and contours** | `distance_to_front()`, `front_frequency()`, `distance_to_contour()`, `distance_to_isobath()`, `distance_to_shore()` |
| **Flow structure** | `eke()`, `current_speed()`, `flow_deformation()`, `ftle()`, `fsle()`, `detect_eddies()`, `distance_to_eddy()`, `residence_time()` |
| **Regional indices** | `scotian_shelf_inflow()`, `northeast_channel_inflow()`, `eastern_gom_salinity()`, `section_transport()`, with `*_section()` and `*_box()` geometries, and `derived_indices()` to list them |
| **Archiving** | `eml_attributes()`, `eml_custom_units()`, `eml_col_classes()` |

Two choices worth knowing before you use any of it. Horizontal gradients come
back in **covariate units per kilometre** rather than per degree, because a
degree of longitude is about 83 km at 42°N and 111 km at the equator, so a
per-degree gradient is stretched by latitude and not comparable across a study
area. And **distance to a front is usually a better predictor than the local
gradient** — a station in a smooth patch a kilometre from a front and one in a
smooth patch a hundred kilometres away have the same gradient and very different
prospects.

→ [**What derivoce computes**](vignettes/derived-covariates.Rmd) is the reference
for every quantity: what it means, why it is computed the way it is, FTLE versus
FSLE, the three ways to measure the same inflow, what the warnings mean, and what
changes when the input has been resampled or gap-filled.

## Column names

Every default is a datamatch catalog name: `SST` and `BOTT` for
`vertical_gradient()`, `UO` and `VO` for `eke()`, `current_speed()`, `ftle()` and
`fsle()`, and `DEPTH` for `distance_to_isobath()`.

Every access function returns columns under the names you asked for, rather than
under the source's own codes, so a catalog fetch needs no column arguments. Here
is the same `eke()` call twice, differing only in the fetch:

```r
bb <- list(xmin = -76, xmax = -65, ymin = 35, ymax = 45)

# Catalog names, so the columns are UO and VO - what eke() expects by default.
env <- datamatch::accessCopernicus(vars = c("UO", "VO"), years = 2010, months = 1:12,
                               bounding_box = bb)
env <- eke(env)

# Copernicus codes, so the columns are uo and vo and have to be named.
env <- datamatch::accessCopernicus(vars = c("uo", "vo"), years = 2010, months = 1:12,
                               bounding_box = bb)
env <- eke(env, u = "uo", v = "vo")
```

`distance_to_isobath()` is the exception: its `DEPTH` column comes from
`datamatch::attach_bathymetry()` rather than from an access function — though
`accessFVCOM()` serves a `DEPTH` of its own, being a model with its own
bathymetry.

## Requirements on the input

Spatial derivatives are only defined on a grid, so the spatial derivations
require points on a **regular lon/lat lattice**. That is the one place a source
matters, and it is about geometry rather than provenance.

| Input | Spatial derivations | Everything else |
|---|---|---|
| Regular lattice — Copernicus, HYCOM, CCMP, most ERDDAP grids | yes | yes |
| Unstructured mesh — `accessFVCOM()` | no | yes |
| Scattered observations | no | yes |

An unstructured mesh returns one row per mesh node, and those nodes are spaced
irregularly by design: resolution follows the coastline instead of a lattice.
Irregular input is rejected rather than interpolated, because a gradient
computed from interpolated data mostly measures the interpolation. Regrid first
with `datamatch::upscale_grid()` or `downscale_grid()` if you need the spatial
derivations on a mesh.

"Everything else" is most of the package: every temporal derivation — lags,
rolling summaries, integrals, anomalies, decomposition, heatwaves — plus
`potential_density()`, `distance_to_shore()`, `box_anomaly()`, `index_series()`
and the EML helpers. Those match points by coordinate and never need a lattice,
so they work on a mesh unchanged.

Central differences are undefined at the grid edge, so boundary cells come back
`NA`. So do the first *n* steps of a lag and the first step of a temporal
gradient. These are genuine absences, not failures.

## Describing the output for an archive

Derived covariates are the hardest part of a dataset to document, because their
meaning lives in how they were computed rather than in what was measured.
`SST_grad` is degrees per kilometre by central differences on a lon/lat lattice,
and nothing in the column name or the numbers says so. That knowledge is already
in this package.

`eml_attributes()` emits it as the attribute table [Ecological Metadata
Language](https://eml.ecoinformatics.org/) wants, ready for
`EML::set_attributes()`. derivoce does not depend on the `EML` package and does
not write XML — it hands over the table:

```r
attributes <- eml_attributes(env)
custom <- eml_custom_units(attributes)
```

A derived unit depends on the source unit — a gradient of temperature is °C/km,
a gradient of chlorophyll mg/m³/km — so the source has to be known first.
Everything datamatch serves already is: the Copernicus physics and
biogeochemistry variables, the seafloor terrain from `attach_bathymetry()`, and
the climate indices from `attach_climate_index()`. A workflow built on those
needs no `units` argument. Pass one for anything else, or to override a default.

Anything still unresolved comes back as `NA` rather than a guess: an archived
dataset with confidently wrong units is worse than one with an obvious hole.
`eml_custom_units()` then returns the declarations EML's fixed 195-entry
dictionary lacks, so the document validates.

## Still to come

- **Vertical gradients from a full depth profile.** The two-level case is
  covered — `vertical_gradient()` takes any two temperature columns and
  `buoyancy_frequency()` gives N² between any two depths. A true profile, and the
  potential energy anomaly that needs one, is still assembly work:
  `accessCopernicus()` returns one level per call, so a profile is several fetches
  joined as columns. Doing that inside the functions is what remains.
- **Extending the LCR index past 2014.** Tried twice and shelved. Monthly
  Copernicus fields fail for a physical reason: the averaging removes the narrow
  Labrador Current jet and the Grand Banks bifurcation the index depends on.
  Daily fields with OceanParcels clear that obstacle and still do not reproduce
  the published series — and it is not the arrival regions, the particle count or
  the domain, all of which were tested and eliminated. What is left is that this
  counts particles entering a region while the paper counts them crossing a
  hydrographic section, and those coordinates are not published.
  [`docs/lcr-extension-experiment.md`](docs/lcr-extension-experiment.md) records
  both attempts.
- **Gulf Stream Index.** NAO, AO, AMO, PDO, LCR, and AMOC are all in
  [datamatch](https://github.com/chross22/datamatch) via
  `attach_climate_index()`. The Gulf Stream Index is harder: it has several
  competing definitions published in papers rather than at a stable URL, so it
  needs a decision about which one.

## Documentation

| | |
|---|---|
| [Getting started](vignettes/derivoce.Rmd) | a worked example end to end, plus the full reference list |
| [What derivoce computes](vignettes/derived-covariates.Rmd) | every quantity, what it means, and why it is computed that way |
| [`docs/methods.md`](docs/methods.md) | longer form, with the reasoning behind each quantity |

## References

Every derived quantity that follows a published method carries that method's
citation. The full list — data sources, software, and the papers behind each
covariate — is at the end of the [getting started
vignette](vignettes/derivoce.Rmd). Each function's own `?help` carries the
references relevant to it, and `as.data.frame(derived_indices())$source` names
the source of every index at runtime.

A scheduled workflow re-checks these each quarter: that every DOI is still
registered, that everything cited in the code or docs appears in the list, and
that nothing in the list is cited nowhere. It opens an issue rather than trying
to fix anything, since choosing the right replacement reference is a judgement
rather than a lookup.

### Citing derivoce

`citation("derivoce")` gives the current form. If a specific covariate follows a
published method, cite that paper too: the list says which, and each function's
`?help` repeats it.
