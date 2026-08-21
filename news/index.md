# Changelog

## derivoce (development version)

### Works with every datamatch source

datamatch now serves seven sources — Copernicus Marine, HYCOM, CCMP
winds, FVCOM, ERDDAP, seafloor terrain and climate indices — and the
documentation here described only the first. It now describes the shape
rather than the product, which is what the package has always actually
required: one row per location and time step.

- The one real distinction is grid geometry, not provenance. Spatial
  derivations need a regular lon/lat lattice, which Copernicus, HYCOM,
  CCMP and most ERDDAP grids give. `accessFVCOM()` returns an
  unstructured mesh, whose nodes are irregularly spaced by design, so
  those derivations refuse it. Every temporal derivation, plus
  [`potential_density()`](https://camilleross.org/derivoce/reference/potential_density.md),
  [`distance_to_shore()`](https://camilleross.org/derivoce/reference/distance_to_shore.md),
  [`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md),
  [`index_series()`](https://camilleross.org/derivoce/reference/index_series.md)
  and the EML helpers, works on a mesh unchanged. README and
  [`?rasterize_step`](https://camilleross.org/derivoce/reference/rasterize_step.md)
  set this out, and tests pin it.
- That refusal used to arrive as terra’s “error in evaluating the
  argument ‘x’”, because the grid was passed to the derivation as a
  promise and failed inside terra rather than before it. It is now
  forced first, so the explanation survives — and the explanation now
  names the mesh case and points at
  [`datamatch::upscale_grid()`](https://camilleross.org/datamatch/reference/upscale_grid.html).
- [`eml_attributes()`](https://camilleross.org/derivoce/reference/eml_attributes.md)
  knows the units of all 42 variables across the seven sources,
  including CCMP winds, HYCOM bottom velocities, and FVCOM’s
  depth-averaged velocities, surface fluxes and wind stress.
- Two things this documentation had recorded as impossible are not any
  more: wind stress is served (`TAUX`, `TAUY`, `TAU`), and so is bottom
  salinity (`BOTS`), which makes a surface-to-bottom *density*
  difference available where only a temperature difference was before.
  `docs/methods.md` records both corrections rather than quietly editing
  them away.

Seven new derivations, none of which need any input `datamatch`’s access
functions do not already serve.

### Anomalies and extremes

- [`cell_anomaly()`](https://camilleross.org/derivoce/reference/cell_anomaly.md)
  removes each cell’s own mean, leaving the departure rather than the
  geography — the cell-wise counterpart of
  [`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md).
  `standardize = TRUE` gives a z-score, comparable across places of
  different variability at the cost of the magnitude.
- [`cell_anomaly()`](https://camilleross.org/derivoce/reference/cell_anomaly.md)
  gains `detrend`, which removes a fitted linear trend as well as the
  mean. In a warming shelf sea an anomaly that still contains the trend
  largely encodes which year it is, and a model given it will fit the
  trend and appear to have learned something about temperature.
- [`decompose_covariate()`](https://camilleross.org/derivoce/reference/decompose_covariate.md)
  splits each cell’s series into a centred trend, a repeating seasonal
  cycle, a residual, and optionally the slope in units per year,
  additively so that `value = mean + trend + seasonal + residual`. The
  trend and the cycle are estimated in one fit rather than sequentially:
  taken one after the other, each absorbs part of the other, and a
  series containing nothing but a seasonal cycle acquires a spurious
  trend worth a few percent of its amplitude.
- [`index_series()`](https://camilleross.org/derivoce/reference/index_series.md)
  collapses a broadcast per-step column back to one row per time step,
  for plotting or export. It refuses a column that varies within a step,
  since collapsing a map would keep one arbitrary cell and discard the
  pattern.
- [`marine_heatwave()`](https://camilleross.org/derivoce/reference/marine_heatwave.md)
  implements Hobday et al. (2016, 2018): a percentile threshold on each
  cell’s own seasonal climatology, events as runs of consecutive
  exceeding steps, and intensity, duration, cumulative intensity and the
  four categories. `direction = "cold"` gives cold spells.

Both fail silently given too little history — a monthly climatology over
one year is exactly zero everywhere, and a 90th percentile over three
values is the largest of the three — so both warn and say which way out
applies.

### Flow structure

- [`flow_deformation()`](https://camilleross.org/derivoce/reference/flow_deformation.md)
  gives vorticity, divergence, the two strain components and their
  magnitude, the Okubo-Weiss parameter, and the Rossby number. It fills
  the gap between
  [`eke()`](https://camilleross.org/derivoce/reference/eke.md), which
  needs a series, and
  [`ftle()`](https://camilleross.org/derivoce/reference/ftle.md)/[`fsle()`](https://camilleross.org/derivoce/reference/fsle.md),
  which need trajectories.
- [`detect_eddies()`](https://camilleross.org/derivoce/reference/detect_eddies.md)
  identifies eddies as objects rather than as a field, by the
  Okubo-Weiss criterion of Isern-Fontanet et al. (2003), and reports
  whether a cell is inside one, which way it turns, and how big it is.
  Polarity matters ecologically: cyclonic cores upwell and anticyclonic
  ones downwell, so a covariate that only says “eddy” averages two
  opposite things together.
- [`distance_to_eddy()`](https://camilleross.org/derivoce/reference/distance_to_eddy.md)
  joins the `distance_to_*` family, optionally restricted to one
  polarity.
- [`residence_time()`](https://camilleross.org/derivoce/reference/residence_time.md)
  releases a particle at every point in a box and measures how long it
  stays, forward or backward. Right-censored at `max_days`, which the
  documentation and a warning both insist on: averaging the column
  biases it downwards, most severely at the most retentive sites.
- [`front_frequency()`](https://camilleross.org/derivoce/reference/front_frequency.md)
  measures how reliably a place is frontal, where
  [`distance_to_front()`](https://camilleross.org/derivoce/reference/distance_to_front.md)
  measures how far one was at a moment.

### Water properties and history

- [`potential_density()`](https://camilleross.org/derivoce/reference/potential_density.md)
  is the UNESCO (1983) one-atmosphere equation of state, reproducing the
  published check values to within 5e-6 kg/m^3.

- [`buoyancy_frequency()`](https://camilleross.org/derivoce/reference/buoyancy_frequency.md)
  computes N^2 between two depths, the stratification measure
  [`vertical_gradient()`](https://camilleross.org/derivoce/reference/vertical_gradient.md)
  only approximates with a temperature difference. That matters where
  salinity varies: Scotian Shelf inflow is fresh enough to stratify
  water barely warmer at the surface.

- [`eady_growth_rate()`](https://camilleross.org/derivoce/reference/eady_growth_rate.md)
  gives the growth rate of baroclinic instability, after Eady (1949) in
  the maximum-growth form of Lindzen and Farrell (1980). It complements
  [`detect_eddies()`](https://camilleross.org/derivoce/reference/detect_eddies.md)
  — that finds eddies that exist, this finds where conditions favour
  making them. Eady is a person, not a spelling of “eddy”.

  Both need velocities or densities at two levels, which is two
  `datamatch::accessCopernicus()` calls at different `depth` ranges
  joined as columns. A previous version of `docs/methods.md` asserted
  this was impossible because only surface values and a bottom
  temperature were available; that was wrong, and is corrected there.

- [`rolling_covariate()`](https://camilleross.org/derivoce/reference/rolling_covariate.md)
  summarises a trailing window — mean, sd, min, max, sum, median, range
  — with the same step-versus-calendar distinction as
  [`lag_covariate()`](https://camilleross.org/derivoce/reference/lag_covariate.md),
  on the same month counter.

### Metadata

- [`eml_attributes()`](https://camilleross.org/derivoce/reference/eml_attributes.md)
  describes the derived columns as an Ecological Metadata Language
  attribute table, in the shape `EML::set_attributes()` consumes.
  Derived covariates are the hardest part of a dataset to document,
  because their meaning lives in how they were computed rather than in
  what was measured, and that knowledge is already here.
- [`eml_custom_units()`](https://camilleross.org/derivoce/reference/eml_custom_units.md)
  returns declarations for the units EML’s dictionary does not carry —
  per second, per second squared, per day, and metres squared per second
  squared among them. A document using an undeclared unit does not
  validate, so a test asserts that every unit the package can emit is
  either standard or declared.
- [`eml_col_classes()`](https://camilleross.org/derivoce/reference/eml_col_classes.md)
  supplies the matching `col_classes` vector.
- Units for everything `datamatch` serves are built in — the Copernicus
  physics and biogeochemistry variables, the seafloor terrain from
  `attach_bathymetry()`, and the climate indices from
  `attach_climate_index()` — so a workflow built on those needs no
  `units` argument, and a derived unit such as a gradient’s is composed
  from its source. A test checks the table against datamatch’s live
  catalogue whenever datamatch is installed, so it cannot quietly fall
  behind. It has already earned its keep: it caught seven variables
  added to datamatch part way through this work — bottom salinity
  (`BOTS`), 10 m wind (`WSPD`, `UWND`, `VWND`) and wind stress (`TAUX`,
  `TAUY`, `TAU`) — which had made a statement in `docs/methods.md` about
  wind stress being unavailable go stale within the day.

No dependency on the `EML` package: derivoce hands over the table and
lets EML write the XML. Units that depend on the source column, such as
a gradient’s, come back as `NA` unless `units` says what the source
holds — a guessed unit in an archived dataset is worse than a visible
gap.

### Documentation

- `docs/lcr-extension-experiment.md` now records the second attempt at
  extending the Labrador Current retroflection index. Daily fields with
  OceanParcels clear the obstacle that stopped the monthly attempt and
  still do not reproduce the published series; the arrival regions, the
  particle count and the domain were each tested and eliminated as
  explanations.

## derivoce 0.1.0

First release.

derivoce takes the output of
[`datamatch::accessEnvDat()`](https://camilleross.org/datamatch/reference/accessEnvDat.html)
— an `sf` point object per time step — and returns the same shape with
derived columns added, so derived covariates flow into a species
distribution model alongside the variables they came from.

### Gradients and change over time

- [`horizontal_gradient()`](https://camilleross.org/derivoce/reference/horizontal_gradient.md)
  and
  [`vertical_gradient()`](https://camilleross.org/derivoce/reference/vertical_gradient.md)
  for spatial gradients, and
  [`temporal_gradient()`](https://camilleross.org/derivoce/reference/temporal_gradient.md)
  for rate of change between time steps.
- [`integrate_covariate()`](https://camilleross.org/derivoce/reference/integrate_covariate.md)
  accumulates a variable over a time window.
- [`lag_covariate()`](https://camilleross.org/derivoce/reference/lag_covariate.md)
  lags a covariate either by position (`by = "step"`) or by calendar
  time (`by = "day"`, `"month"`, `"year"`). The two agree until the
  record has a gap and then disagree silently, so the distinction is
  explicit rather than assumed: “three months ago” is a statement about
  the organism, “three steps ago” is a statement about the sampling.

### Fluid dynamics

- [`eke()`](https://camilleross.org/derivoce/reference/eke.md) for eddy
  kinetic energy and
  [`current_speed()`](https://camilleross.org/derivoce/reference/current_speed.md)
  for speed from velocity components.
- [`ftle()`](https://camilleross.org/derivoce/reference/ftle.md) and
  [`fsle()`](https://camilleross.org/derivoce/reference/fsle.md) for
  finite-time and finite-size Lyapunov exponents, advecting particles
  with fourth-order Runge-Kutta through a time-varying field.
- Both Lyapunov exponents lose a margin at the domain edge, where the
  neighbouring particles a deformation needs fall outside the data. The
  size of that margin is documented on each function, and empty output
  is explained rather than returned silently.

### Distances

- [`distance_to_shore()`](https://camilleross.org/derivoce/reference/distance_to_shore.md),
  [`distance_to_isobath()`](https://camilleross.org/derivoce/reference/distance_to_isobath.md),
  [`distance_to_front()`](https://camilleross.org/derivoce/reference/distance_to_front.md)
  and
  [`distance_to_contour()`](https://camilleross.org/derivoce/reference/distance_to_contour.md).

### Region-scale indices

Indices describing an area rather than a cell:

- [`section_transport()`](https://camilleross.org/derivoce/reference/section_transport.md)
  for volume transport across a section.
- [`water_mass_fraction()`](https://camilleross.org/derivoce/reference/water_mass_fraction.md)
  for water-mass fractions from temperature-salinity endmember mixing.
- [`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md)
  for regional property anomalies.
- Named cases for the regions this was built for:
  [`scotian_shelf_inflow()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow.md),
  [`northeast_channel_inflow()`](https://camilleross.org/derivoce/reference/northeast_channel_inflow.md)
  and
  [`eastern_gom_salinity()`](https://camilleross.org/derivoce/reference/eastern_gom_salinity.md),
  with their geometry exposed as
  [`scotian_shelf_inflow_section()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow_section.md),
  [`northeast_channel_section()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow_section.md)
  and
  [`eastern_gom_box()`](https://camilleross.org/derivoce/reference/eastern_gom_box.md)
  so the definitions can be inspected and plotted.
- [`derived_indices()`](https://camilleross.org/derivoce/reference/derived_indices.md)
  catalogues them: what each measures, what it needs as input, and where
  the concept comes from. Scotian Shelf inflow appears three times
  because the literature measures it three ways — a transport, a
  water-mass fraction and a box anomaly — and they answer different
  questions.

The named sections were placed against real Copernicus velocities rather
than by eye, and the Northeast Channel section on 60 months rather than
a single season. How that was done is recorded in `docs/`.

### Documentation

- A vignette that runs its own examples, so what it shows is what the
  code currently does.
- `docs/methods.md` for the derivations behind the functions.
- Citations for the derived covariates and the source products they come
  from, on each function’s help page and in the README. `inst/CITATION`
  takes its version from `DESCRIPTION` rather than hard-coding one, and
  a test enforces that.
- A quarterly check that the cited links have not rotted.

### Known limitations

- Covariates that cannot be informed by the data given — too few time
  steps for a lag, too small a domain for a Lyapunov exponent — warn
  rather than returning a column of `NA` without comment.
- [`datamatch::attach_climate_index()`](https://camilleross.org/datamatch/reference/attach_climate_index.html)
  serves `LCR`, the Labrador Current retroflection index, only for
  1993–2014. Two attempts to recompute it so the series could be
  extended are recorded in `docs/lcr-extension-experiment.md`; neither
  reproduces the published index, and the reasons are documented there.
