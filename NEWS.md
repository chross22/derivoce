# derivoce (development version)

Seven new derivations, none of which need any input
`datamatch::accessEnvDat()` does not already serve.

## Anomalies and extremes

* `cell_anomaly()` removes each cell's own mean, leaving the departure rather
  than the geography — the cell-wise counterpart of `box_anomaly()`.
  `standardize = TRUE` gives a z-score, comparable across places of different
  variability at the cost of the magnitude.
* `cell_anomaly()` gains `detrend`, which removes a fitted linear trend as well
  as the mean. In a warming shelf sea an anomaly that still contains the trend
  largely encodes which year it is, and a model given it will fit the trend and
  appear to have learned something about temperature.
* `decompose_covariate()` splits each cell's series into a centred trend, a
  repeating seasonal cycle, a residual, and optionally the slope in units per
  year, additively so that
  `value = mean + trend + seasonal + residual`. The trend and the cycle are
  estimated in one fit rather than sequentially: taken one after the other, each
  absorbs part of the other, and a series containing nothing but a seasonal
  cycle acquires a spurious trend worth a few percent of its amplitude.
* `index_series()` collapses a broadcast per-step column back to one row per
  time step, for plotting or export. It refuses a column that varies within a
  step, since collapsing a map would keep one arbitrary cell and discard the
  pattern.
* `marine_heatwave()` implements Hobday et al. (2016, 2018): a percentile
  threshold on each cell's own seasonal climatology, events as runs of
  consecutive exceeding steps, and intensity, duration, cumulative intensity and
  the four categories. `direction = "cold"` gives cold spells.

Both fail silently given too little history — a monthly climatology over one
year is exactly zero everywhere, and a 90th percentile over three values is the
largest of the three — so both warn and say which way out applies.

## Flow structure

* `flow_deformation()` gives vorticity, divergence, the two strain components
  and their magnitude, the Okubo-Weiss parameter, and the Rossby number. It
  fills the gap between `eke()`, which needs a series, and `ftle()`/`fsle()`,
  which need trajectories.
* `residence_time()` releases a particle at every point in a box and measures
  how long it stays, forward or backward. Right-censored at `max_days`, which
  the documentation and a warning both insist on: averaging the column biases
  it downwards, most severely at the most retentive sites.
* `front_frequency()` measures how reliably a place is frontal, where
  `distance_to_front()` measures how far one was at a moment.

## Water properties and history

* `potential_density()` is the UNESCO (1983) one-atmosphere equation of state,
  reproducing the published check values to within 5e-6 kg/m^3.
* `rolling_covariate()` summarises a trailing window — mean, sd, min, max, sum,
  median, range — with the same step-versus-calendar distinction as
  `lag_covariate()`, on the same month counter.

## Documentation

* `docs/lcr-extension-experiment.md` now records the second attempt at
  extending the Labrador Current retroflection index. Daily fields with
  OceanParcels clear the obstacle that stopped the monthly attempt and still do
  not reproduce the published series; the arrival regions, the particle count
  and the domain were each tested and eliminated as explanations.

# derivoce 0.1.0

First release.

derivoce takes the output of `datamatch::accessEnvDat()` — an `sf` point object
per time step — and returns the same shape with derived columns added, so
derived covariates flow into a species distribution model alongside the
variables they came from.

## Gradients and change over time

* `horizontal_gradient()` and `vertical_gradient()` for spatial gradients, and
  `temporal_gradient()` for rate of change between time steps.
* `integrate_covariate()` accumulates a variable over a time window.
* `lag_covariate()` lags a covariate either by position (`by = "step"`) or by
  calendar time (`by = "day"`, `"month"`, `"year"`). The two agree until the
  record has a gap and then disagree silently, so the distinction is explicit
  rather than assumed: "three months ago" is a statement about the organism,
  "three steps ago" is a statement about the sampling.

## Fluid dynamics

* `eke()` for eddy kinetic energy and `current_speed()` for speed from velocity
  components.
* `ftle()` and `fsle()` for finite-time and finite-size Lyapunov exponents,
  advecting particles with fourth-order Runge-Kutta through a time-varying
  field.
* Both Lyapunov exponents lose a margin at the domain edge, where the
  neighbouring particles a deformation needs fall outside the data. The size of
  that margin is documented on each function, and empty output is explained
  rather than returned silently.

## Distances

* `distance_to_shore()`, `distance_to_isobath()`, `distance_to_front()` and
  `distance_to_contour()`.

## Region-scale indices

Indices describing an area rather than a cell:

* `section_transport()` for volume transport across a section.
* `water_mass_fraction()` for water-mass fractions from temperature-salinity
  endmember mixing.
* `box_anomaly()` for regional property anomalies.
* Named cases for the regions this was built for:
  `scotian_shelf_inflow()`, `northeast_channel_inflow()` and
  `eastern_gom_salinity()`, with their geometry exposed as
  `scotian_shelf_inflow_section()`, `northeast_channel_section()` and
  `eastern_gom_box()` so the definitions can be inspected and plotted.
* `derived_indices()` catalogues them: what each measures, what it needs as
  input, and where the concept comes from. Scotian Shelf inflow appears three
  times because the literature measures it three ways — a transport, a
  water-mass fraction and a box anomaly — and they answer different questions.

The named sections were placed against real Copernicus velocities rather than
by eye, and the Northeast Channel section on 60 months rather than a single
season. How that was done is recorded in `docs/`.

## Documentation

* A vignette that runs its own examples, so what it shows is what the code
  currently does.
* `docs/methods.md` for the derivations behind the functions.
* Citations for the derived covariates and the source products they come from,
  on each function's help page and in the README. `inst/CITATION` takes its
  version from `DESCRIPTION` rather than hard-coding one, and a test enforces
  that.
* A quarterly check that the cited links have not rotted.

## Known limitations

* Covariates that cannot be informed by the data given — too few time steps for
  a lag, too small a domain for a Lyapunov exponent — warn rather than
  returning a column of `NA` without comment.
* `datamatch::attach_climate_index()` serves `LCR`, the Labrador Current
  retroflection index, only for 1993–2014. Two attempts to recompute it so the
  series could be extended are recorded in
  `docs/lcr-extension-experiment.md`; neither reproduces the published index,
  and the reasons are documented there.
