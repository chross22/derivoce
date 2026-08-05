# How each derived covariate is computed, and why

This is the reasoning behind the numbers: what each quantity means, how it is
calculated, which choices were deliberate, and where the results differ from the
older pipeline they replace. The function documentation says *what*; this says
*why*.

## Context

The pipeline these covariates feed —
[`taupatch`](https://github.com/chross22/taupatch), rebuilt from Ross et al.
(2023) — models high-abundance zooplankton patches from environmental
conditions. The environmental variables Copernicus serves directly (temperature,
salinity, currents, chlorophyll) describe *what conditions are*. The derived
quantities here describe *how conditions are arranged and how they are changing*,
which is often what actually concentrates plankton:

- A **front** — a sharp horizontal gradient — is a convergence zone where
  buoyant particles and weak swimmers accumulate. The temperature itself may be
  unremarkable on both sides; the *change* is what matters.
- **Stratification** — a vertical gradient — determines whether the surface layer
  is isolated from deep water, which controls both nutrient supply and where
  animals can stay.
- **Accumulated** food matters more than instantaneous food, because a copepod
  sampled in June grew on the spring bloom.
- **Lagged** conditions matter because population responses take time.

## The data shape

Everything takes and returns the output shape of
[`datamatch::accessEnvDat()`](https://github.com/chross22/datamatch): an `sf`
POINT object with one row per (grid point, time step), a column per covariate,
and `YEAR`/`MONTH`/`DAY` columns.

Keeping input and output shapes identical is what lets the functions compose:

```r
env |> horizontal_gradient("thetao") |> lag_covariate("thetao")
```

Internally, spatial operations need a grid, so each time step is rasterized,
operated on, and sampled back to the points (`per_time_step()` in `R/grid.R`).
Temporal operations work directly on the point table, matching locations by
coordinate.

---

## Horizontal gradients

**What:** the magnitude of a covariate's rate of change with horizontal distance,
`|∇C|`, in covariate units per kilometre.

**How:** central differences on the covariate's own grid.

```
∂C/∂x ≈ (C[i+1,j] - C[i-1,j]) / 2Δx
∂C/∂y ≈ (C[i,j+1] - C[i,j-1]) / 2Δy
|∇C|  = sqrt((∂C/∂x)² + (∂C/∂y)²)
```

Central differences rather than one-sided ones because they are second-order
accurate — the leading error term cancels — for the same computational cost. The
price is that the grid boundary has no value, since one neighbour is missing.
Those cells come back `NA` rather than being filled by a one-sided estimate,
because a boundary value computed by a different formula is not comparable to the
interior ones.

### The unit question, which is the whole design

The naive implementation leaves Δx and Δy in **degrees**. That is wrong in a way
that is easy to miss, because the output still looks like a sensible field.

A degree of latitude is about 111.32 km everywhere. A degree of *longitude* is
111.32 km at the equator and shrinks as `cos(latitude)`:

| Latitude | km per degree longitude |
|---|---|
| 0° | 111.3 |
| 35° | 91.2 |
| 42° | 82.7 |
| 45° | 78.7 |

Across taupatch's study area (35°N–45°N) that is a 16% spread. A gradient left in
per-degree units is therefore stretched by latitude: the same physical front
registers as a *weaker* east-west gradient in the north than in the south, purely
as an artifact. A model would learn that artifact as if it were a latitudinal
gradient in front strength.

So the spacing is converted to distance before dividing:

```
Δx(metres) = Δx(degrees) × 111320 × cos(latitude)
Δy(metres) = Δy(degrees) × 111320
```

`cos(latitude)` is evaluated **per grid row**, not once for the whole grid. Using
a single central latitude would leave a residual error across the domain —
roughly 30% at 45° for a wide grid.

`cell_size()` handles projected rasters too: if the CRS is already in linear
units, the resolution is used directly with no conversion.

The 111320 m figure is the mean-Earth-radius value. The true value varies by a
few tenths of a percent with latitude because the Earth is oblate; that is far
below the error in the underlying model fields and is not worth correcting.

### Why this is not `raster::terrain()`

The older pipeline computed `sst_grad` and `uv_grad` with `raster::terrain(sst)`.
That function is built for **elevation** models. It assumes the cell values are a
length in the same units as the coordinates and returns a **slope angle**:

```
slope = atan(sqrt((∂z/∂x)² + (∂z/∂y)²))
```

For a temperature field this is dimensionally meaningless — it takes the
arctangent of a quantity in °C per degree-of-longitude, which is not an angle and
not anything else. It also saturates: `atan` compresses everything above a
moderate gradient toward π/2, so genuinely sharp fronts become indistinguishable
from moderate ones.

In practice the old values were probably still *useful*, because `atan` is
monotonic — a bigger gradient still gave a bigger number, so the random forest
could use it as a front detector. But the values were not a rate of change, were
not comparable across latitudes, and could not be interpreted or compared to
anything published in physical units.

**This means derived gradients here are not numerically comparable to the ones in
Ross et al. (2023).** Lags and integrals are; gradients are not. If bit-level
reproduction of the published model is ever needed, that would require adding a
`method = "terrain"` option rather than reinterpreting the current output.

---

## Vertical gradients

**What:** surface temperature minus bottom temperature in each cell — a
stratification index. Large where a warm surface layer sits over cold deep water;
near zero where tide and wind have mixed the column through.

**How:** a column-wise difference, `SST - BOTT`.

There is no depth-profile interpolation here, and deliberately so. Copernicus
serves `thetao` (surface) and `bottomT` (sea floor) as two separate variables in
the *same* physics dataset, so both arrive in a single fetch and the difference is
already the quantity of interest.

This is worth stating because the obvious alternative is a genuine trap:
`datamatch::accessEnvDat()` assigns its output column names **positionally**
(`c("x", "y", vars, "YEAR", "MONTH", "DAY")`). If a depth *range* were requested
spanning several model levels, the raster would come back with more layers than
there are variable names, and that assignment would either error or silently
mislabel columns. Requesting a multi-level profile to compute `dT/dz` properly
would need a datamatch change first. The surface-to-bottom difference sidesteps
that entirely.

With a `depth` column supplied, the difference is divided by depth to give a
per-metre rate rather than a total difference. Depths of zero or less become `NA`
rather than `Inf` or a sign flip — a cell on land or at the waterline has no
meaningful stratification rate.

---

## Temporal gradients

**What:** the rate of change at a fixed location between consecutive time steps.
A water mass warming quickly is a different habitat from one sitting at the same
temperature.

**How:** `(C[t] - C[t-1]) / Δt`, built on `lag_covariate()`.

`Δt` depends on `per`:

- `"step"` (default) — change per time step, `Δt = 1`. Correct when steps are
  evenly spaced, which monthly products are in index terms.
- `"day"` — divides by the actual number of days between steps, from the calendar
  dates. Month lengths vary from 28 to 31 days, an 11% spread, so this differs
  from `"step"` by more than rounding.
- `"month"` — divides by days, then by 30.4375 (the mean Gregorian month), giving
  a nominal per-month rate on an even footing.

The first time step has no predecessor and is `NA`.

---

## Lags

**What:** a covariate's value `n` time steps earlier at the same location.

**How:** for each step, find the step `n` back and match points by coordinate.

**Matching by coordinate, not by row order,** is the important detail. It would
be simpler to assume both steps list their grid points in the same sequence, and
that assumption usually holds — until a fetch is reordered, a step is missing
points, or two products are joined. When it breaks, the failure is silent: every
row gets *a* value, just the wrong cell's. So `location_key()` builds a rounded
coordinate string and matches on that. Rounding (6 decimal places, ~0.1 m)
absorbs float noise from coordinates written and re-read at different precisions.

The first `n` steps have no predecessor and are `NA`.

`n = 1` reproduces the older pipeline's `lag_sst` exactly.

---

## Time integration

**What:** a covariate accumulated over preceding time steps at each location. The
food a survey encounters is the food built up since the season began, not the
instantaneous concentration.

**How:** sum over the contributing steps, with `window` controlling which those
are:

- `"year"` (default) — every step so far in the same calendar year, reset at the
  year boundary. This reproduces the older pipeline's `int_chl`, which summed
  chlorophyll from January.
- `"all"` — every step from the start of the record, never resetting.
- a positive integer — a rolling total over that many trailing steps, crossing
  year boundaries.

The year reset is not arbitrary. It encodes the assumption that the biological
year restarts — that January's productivity is the beginning of a new growing
season rather than a continuation of December's. That is defensible in the Gulf
of Maine, where the spring bloom is the dominant annual event, and it is what the
published model did. A rolling window makes no such assumption, which is why it
is available.

**Missing values are counted, not treated as zero.** If a location is absent from
a contributing step, adding zero would understate the total while looking like a
real number. Instead, each contributing step is tracked; a location with no
contributing steps at all returns `NA`, rather than a spurious total of 0.

---

## Distance to fronts and contours

**What:** how far each point is from the nearest front, or from a named contour
such as an isobath.

**Why not just use the gradient:** the local gradient is a blunt predictor of
frontal influence. A station in the middle of a smooth patch has a gradient of
zero whether the nearest front is 2 km away or 200 km away, and those are very
different places to be — plankton accumulate *near* fronts, not only exactly on
them. Distance separates the two cases; the gradient cannot.

**How:** threshold the gradient to identify frontal cells, then run a distance
transform (`terra::distance()`) from every cell to the nearest frontal one.

### What counts as a front

`scope` decides what the quantile threshold is taken over, and the two answers
answer different questions:

- **`"record"`** (default) — one threshold across the whole series, so "front"
  means the same physical sharpness every month. A weakly stratified winter may
  then contain *no* fronts. That is a real result, and suppressing it would be
  the error.
- **`"step"`** — a separate threshold per time step, so every month has fronts by
  construction. Right when the question is where this month's sharpest features
  are; wrong when the question is whether this month has sharp features at all.

A time step with no front gives `NA`, not `0`. Reporting zero would assert "you
are standing on a front", the opposite of the truth.

### Contours and isobaths

`distance_to_contour()` measures to where a covariate crosses a value —
`distance_to_isobath()` is the depth case. Position relative to a contour is
frequently the better predictor: plankton track the shelf break, and "20 km
inshore of the 100 m isobath" locates that far better than "depth = 85 m".

A contour almost never falls exactly on a cell centre, so a cell is marked as on
the contour when the level lies between its value and any neighbour's — detected
with a 3×3 focal range of the above/below indicator. Testing for exact equality
would find nothing on coarsely quantised data.

---

## FTLE

**What:** the Finite-Time Lyapunov Exponent — how fast two initially neighbouring
water parcels separate, in units of 1/day. Ridges mark Lagrangian coherent
structures: the lines that organise surface transport.

**Direction is the substantive choice**, which is why it is an argument:

- **Backward** (default) integrates into the past, so ridges are *attracting*
  structures — where water arriving from different origins converges. This is
  where buoyant particles and weak swimmers pile up, so it is usually the
  direction of interest for aggregation.
- **Forward** integrates into the future, so ridges are *repelling* structures —
  barriers that neighbouring parcels straddle and are pulled apart by.

**How:**

1. Seed particles on the covariate grid at each time step.
2. Advect them with fourth-order Runge-Kutta through the time-varying velocity
   field, bilinear in space and linear in time between bracketing fields.
3. Differentiate the resulting flow map by central differences to get the
   deformation gradient `F`.
4. Form the Cauchy-Green strain tensor `C = FᵀF` and take its largest eigenvalue.
5. `σ = ln(√λ_max) / |T|`.

The 2×2 eigenvalue is closed-form, so there is no per-cell `eigen()` call.

### Details that matter

**Velocities are converted from m/s to degrees per day at every evaluation
point**, with the longitude conversion scaled by `cos(latitude)` — the same
meridian-convergence issue as for gradients, but here it affects the trajectory
itself, not just the reported unit.

This has a consequence worth stating, because it looks like a bug: **a constant
velocity in m/s does not give zero FTLE.** Meridians converge poleward, so a
parcel at 43°N sweeps through more longitude per second than one at 42°N. That
differential is genuine shear. Across 42–43°N it is about 1.6%, giving an FTLE
around 6×10⁻⁴/day — orders of magnitude below real fronts (0.01–0.1/day), but not
zero. Only a motionless field gives exactly zero.

**Positions are converted to metres on an equirectangular plane about the grid
centre** before differentiating, so the strain tensor is computed in one
consistent metric rather than mixing degrees of longitude and latitude, which have
different physical lengths.

**λ < 1 is clamped to 1.** A window over which parcels net *converge* is a real
result, but it is not an exponential separation rate, and `ln` of it would be
negative. Clamping reports zero separation rather than a negative exponent.

**Particles are clamped at the ends of the record** rather than extrapolated: one
that runs past the available fields is held in the last known flow, a milder error
than inventing velocities.

### On monthly data

Monthly-mean velocity fields have already averaged away the eddies that generate
the sharpest structures. FTLE from monthly means is a real diagnostic of the mean
circulation, but it is not the same quantity as FTLE from daily fields, and its
ridges are correspondingly smoother. Prefer a `P1D` product where the choice
exists.

`integration_days` sets the scale of structure resolved: too short and no
separation has accumulated; too long and ridges blur as particles wander far from
where they started. Days to a few weeks is the usual mesoscale range.

---

## Eddy kinetic energy

**What:** `EKE = ½(u′² + v′²)`, the kinetic energy in the *departure* of the flow
from its mean. High EKE marks an energetic, variable region — meanders, rings,
eddies — as distinct from a strong but steady current, which carries high total
kinetic energy and little eddy energy.

**The reference is the whole design.** EKE is only defined relative to a mean, and
that choice decides what counts as an eddy. It is a scientific decision, so it is
an argument:

- **`"record"`** (default) — the mean over the whole series at each location. The
  textbook definition. The seasonal cycle of the mean circulation remains in the
  anomaly, so a region whose currents merely strengthen every summer registers as
  energetic.
- **`"climatology"`** — a separate mean per calendar month, removing the
  repeatable seasonal cycle so only departures from *the usual conditions for that
  month* count. Use when asking about anomalous years rather than about which
  places are energetic.
- **an integer** — a centred rolling mean over that many steps, removing slow
  drift so only variability faster than the window survives.

Means are computed **per location**, so a spatially varying mean circulation is
removed correctly rather than smeared into the anomaly. A domain-wide mean would
leave a large spurious anomaly wherever the mean flow differs from the domain
average — which, on a shelf with opposing currents, is everywhere.

`current_speed()` gives `sqrt(u² + v²)`, reproducing the older pipeline's `uv`;
passing it through `horizontal_gradient()` reproduces `uv_grad`.

---

## Requirements on the input, and what is rejected

`horizontal_gradient()` requires points on a **regular lon/lat lattice**.
`rasterize_step()` verifies this by checking that the unique longitudes and the
unique latitudes are each evenly spaced, and errors if not.

This rejection is deliberate. Gridded model output (Copernicus and similar) is
regular; scattered observations are not. It would be easy to accept anything and
interpolate onto a grid first, and the result would look completely plausible —
but a gradient field computed from interpolated data mostly measures the
interpolation, not the ocean. Failing loudly is better than returning a covariate
that silently encodes the sampling pattern.

Boundary cells are `NA`, as are the first `n` steps of a lag and the first step of
a temporal gradient. These are genuine absences, not failures.

---

## How this is tested

Gradients are checked against fields whose derivatives are known analytically,
which is stronger than checking that output looks reasonable:

- A **uniform field** must give exactly zero everywhere.
- A **north-south ramp** of 1 °C per degree of latitude must give exactly
  `1/111.32` °C/km at *every* point, since latitude spacing does not vary. It must
  also put the entire signal in the northward component and none in the eastward
  one.
- An **east-west ramp** of 1 °C per degree of longitude must give
  `1/(111.32·cos(lat))` — that is, a value that *increases* toward the pole. This
  is the test that catches the per-degree mistake: the naive implementation
  returns a constant here, which is exactly what a plausible-looking wrong answer
  would look like.
- Changing the unit from km to m must scale the result by exactly 1000.
- A field that is flat in space but jumps between months must give zero spatial
  gradient — catching any pooling of time steps.

Temporal functions are checked against fields with a planted linear trend, and
lag matching is checked by **shuffling one time step's rows** and confirming the
answer is unchanged.

---

## Not yet implemented

- **FSLE** — the finite-*size* counterpart to FTLE, integrating until parcels
  reach a chosen separation rather than for a fixed time. It resolves structures
  at a specified spatial scale rather than a specified timescale, which suits
  comparing across regions with different flow speeds. The advection machinery in
  `R/ftle.R` is reusable; the stopping criterion is what differs.
- **Distance to shore** — static, like bathymetry. Needs a coastline source
  (`rnaturalearth` or similar) and a decision about whether "shore" means the
  mainland coast or includes islands.
- **Vertical gradients from a depth profile** — the current implementation is
  surface-minus-bottom. A true `dT/dz` from intermediate levels needs multi-level
  data, which `datamatch::accessEnvDat()` cannot yet return (see above).
- **Climate indices** (NAO, AMO, Gulf Stream Index) — time-only covariates with no
  spatial dimension, broadcast across every point in a time step. These are a
  retrieval problem rather than a derivation one, so they likely belong in
  `datamatch` alongside `accessEnvDat()`.

## FTLE at depth

FTLE is computed from whatever velocity field it is handed, so depth-resolved FTLE
needs only depth-resolved velocities — no change here.

Copernicus does provide them: `GLOBAL_MULTIYEAR_PHY_001_030` (GLORYS12V1) carries
`uo` and `vo` on 50 vertical levels, so a fetch at a chosen `depth` returns the
flow at that level and `ftle()` runs on it unmodified.

The constraint is on the fetching side. `accessEnvDat()` names its output columns
positionally, so it returns **one level per call**. Depth-resolved FTLE therefore
means one fetch per level, each producing its own FTLE field — which works today,
just not in a single call.

Worth knowing before doing it: horizontal velocities weaken and decorrelate with
depth, so sub-surface FTLE ridges are generally weaker and less coherent than
surface ones. Surface FTLE is also the case validated against drifters in the
literature; deeper fields are less constrained by observation.
