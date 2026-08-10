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
