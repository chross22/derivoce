# derivoce audit: product sourcing and uncertainty reporting

Audit date: 2026-10-05. Package version: 0.2.0 (commit dbb5d91).

## 1. Scope

This document records an audit of derivoce covering two questions:

1. How each derived product is sourced: the formula or algorithm, the
   literature it follows, and whether the implementation matches that
   literature.
2. Whether the uncertainty of each product is conveyed to the user: through
   warnings, `NA` handling, output columns, documentation, or not at all.

Files in scope: `R/*.R`, `tests/testthat/`, `docs/methods.md`,
`vignettes/`, `README.md`, `NEWS.md`, `inst/CITATION`, the citation workflow,
and `lcr-pipeline/`. Upstream packages (datamatch, copernicusmarine,
OceanParcels) and the contents of the cited papers were not examined.

## 2. Method

The code was split into four groups and each was read in full by a separate
reviewer: (a) gradients and time series, (b) flow structure, (c) density and
regional products, (d) the LCR pipeline, EML output and citations. Each
reviewer applied the same rubric: provenance, fidelity to the cited method,
uncertainty signalling, and defects with file and line. Defects marked
"reproduced" were confirmed by running a snippet against `pkgload::load_all()`.
The five findings marked "confirmed" in section 7 were re-checked
independently by reading the code, and the rolling-window finding was also run.

## 3. Product sourcing

Fidelity: **exact** = matches the cited method; **partial** = matches with
undocumented deviations; **own** = construction of the package author.
Signalling: **good**, **partial**, or **none**.

| Product | Source cited | Fidelity | Signalling |
|---|---|---|---|
| `horizontal_gradient` | none (central differences) | own, textbook form | partial |
| `vertical_gradient` | none | own (temperature-only proxy) | partial |
| `temporal_gradient`, `lag_covariate`, `integrate_covariate`, `rolling_covariate` | Ross et al. 2023 for lag and integral | own | partial; `integrate_covariate` none |
| `decompose_covariate` | none | own (joint `lm` fit) | partial |
| `cell_anomaly`, `box_anomaly` | none | own | partial |
| `marine_heatwave` | Hobday et al. 2016, 2018 | partial (section 5.2) | partial |
| `front_frequency`, `distance_to_front` | Belkin and O'Reilly 2009 (idea only) | partial; quantile threshold | partial |
| `distance_to_contour`, `distance_to_isobath` | none | own | none |
| `distance_to_shore` | Natural Earth | exact | partial |
| `potential_density` | UNESCO 1983 | exact; 15 coefficients verified | good |
| `buoyancy_frequency` | none cited | partial; fixed rho0, surface-referenced density | partial |
| `eady_growth_rate` | Eady 1949; Lindzen and Farrell 1980 | exact formula | partial |
| `eke` | none | partial; total, not geostrophic, anomaly | none |
| `flow_deformation`, Okubo-Weiss | Okubo 1970; Weiss 1991 | exact; curvature terms omitted | partial |
| `detect_eddies` | Isern-Fontanet et al. 2003 | partial; labelling and radius are own | partial |
| `ftle` | Haller 2015 | exact definition; boundary clamp | partial |
| `fsle` | d'Ovidio et al. 2004 | partial; quantised separation time | partial |
| `residence_time` | none | own | partial |
| `water_mass_fraction` | Townsend et al. 2015 (concept) | own; endmembers uncited | partial |
| `section_transport`, `scotian_shelf_inflow`, `northeast_channel_inflow` | Ramp et al. 1985; Feng et al. 2016; Wang et al. 2022; Du et al. 2022; Silver et al. 2023 (context) | own; endpoints tuned on GLORYS 2008-2012 | partial |
| `eastern_gom_salinity` | Grodsky et al. 2025 | partial; box approximate | partial |
| `eml_attributes` | EML 2.2 unit dictionary | partial | none (section 5.4) |

External data enter the package at three points. None records a version in its
output.

| Item | Origin | Recorded in output |
|---|---|---|
| Velocity and scalar fields | datamatch (Copernicus, HYCOM, CCMP, FVCOM, ERDDAP) | no |
| `lcr_published.csv` | exported from datamatch; table of origin not stated | no |
| Shoreline | Natural Earth via `rnaturalearth`, resolution set by argument | no (`shore_dist` carries no resolution tag) |

## 4. Uncertainty reporting: cross-cutting findings

1. **No product returns a spread.** No function emits a standard error,
   confidence interval, sample count, or quality flag column. `slope`,
   `water_mass_fraction`, transport, front frequency and the climatology
   baselines are all point estimates.
2. **Warnings depend on `NA`, and several failure modes are not `NA`.**
   Frozen-flow advection past the record end (section 5.1), zero anomalies
   from one-year climatologies, and look-ahead windows all return finite,
   plausible values and so trigger no warning.
3. **Warning thresholds are coarse.** Lagrangian warnings fire at 90% `NA`;
   instability and Eady warnings at 10% of cells; the climatology warning
   tests the largest group rather than the smallest.
4. **Input provenance is not detected.** Resampled, gap-filled, monthly-mean
   and irregularly spaced inputs are accepted. The caveats exist in
   `docs/methods.md` and the vignette, not in `?help` or at runtime.
5. **Units are assumed.** `u` and `v` are assumed to be in m/s, temperature in
   degrees Celsius, depth positive down. None is checked.
6. **The LCR pipeline reports its failure well.** The README, `STATUS.md` and
   `lcr-extension-experiment.md` state that the recomputed series does not
   reproduce the published one, and `validate.py` withholds a verdict on short
   records. Exceptions are listed in section 5.5.

## 5. Defects

### 5.1 High severity

| # | Location | Defect | Correction |
|---|---|---|---|
| 1 | `R/transport.R:99-103` | Dropped (land or `NA`) sample points enter the sum as zero. A uniform 1 m/s flow with half the section masked returned 0.58 of the full transport. The docstring at lines 30-33 states the opposite. With the default `min_coverage = 0.5` the bias reaches 50% silently. (reproduced; confirmed) | Rescale by `1 / mean(usable)` or warn below full coverage; add a partial-coverage test. |
| 2 | `R/ftle.R:217`; reached from `fsle.R`, `residence.R` | Time is clamped to the record, so a particle past either end is advected in a frozen flow and returns a finite value. A 3-step record with `integration_days = 14` gave 78% non-`NA` cells and no warning. (reproduced; confirmed) | Return `NA` or warn with the number of affected steps; document in roxygen, not only `methods.md:760`. |
| 3 | `R/rolling.R:155-157` | Calendar windows ignore `DAY`. On daily data every later day of the current month is in the window; `n = 1, by = "month", stat = "max"` returned 31 on 1 January. This puts future data in a covariate. (reproduced; confirmed) | Build the window from timestamps as `by = "day"` does. |
| 4 | `R/region.R:65-71` | `box_anomaly(reference = "climatology")` on a one-year record yields an exact zero anomaly. (reproduced; confirmed by code) | Error or warn when any calendar month has one observation. |
| 5 | `R/eml-registry.R:21-29` | Transport, Scotian inflow and channel inflow are archived as `cubicMetersPerSecond`; the quantity is m^2/s per unit depth (`transport.R:3-4`). (confirmed) | Declare a custom unit, `metersSquaredPerSecond`. |
| 6 | `R/eml-registry.R:149-151` | `_tgrad` is archived with the source unit; it is a rate per step, day or month. (confirmed) | Append the `per` denominator. |
| 7 | `R/eml.R:139-145`, `R/eml-registry.R` | Composed units (`UO_grad`, `CHL_grad`, `CHL_slope`) are not in the registry, and `eml_custom_units()` omits them, so the EML document fails validation. The README's own chlorophyll example is affected. | Emit declarations for every composed unit. |
| 8 | `R/eml.R:179-184` | Fixed-name registry columns (`speed`, `N`, `radius`, `density`) are matched before `units` is consulted. A user column of that name receives a confident unit. This contradicts the README claim that unresolved units return `NA`. (confirmed) | Consult `units` first, or require a derivoce marker. |
| 9 | `R/transport.R:295-321`; `docs/section-placement-diagnostics.R` | Section endpoints were selected in-sample on 2008-2012 by requiring positive net inflow; a rotation of the Cape Sable line changes the index about nine-fold. No out-of-sample test or sensitivity band exists. | State the selection in the function docs; add an endpoint-jitter sensitivity result. |

### 5.2 Medium severity

| # | Location | Defect | Correction |
|---|---|---|---|
| 10 | `R/heatwave.R:181-214`, `29-31` | `max_gap` bridging is applied before `min_steps`; Hobday applies the minimum duration first. Flags `T T F T T F F F` with `min_steps = 5, max_gap = 2` give a 5-step event. The claim that these settings reproduce the published definition is false. (reproduced) | Reorder operations or reword the claim. |
| 11 | `R/heatwave.R` | Climatology and threshold are monthly step functions, not an 11-day windowed day-of-year climatology over a 30-year baseline. Category is assigned per row, not per event peak. No fixed-baseline argument exists although the docs describe one. | Document each deviation in `?marine_heatwave`. |
| 12 | `R/anomaly.R:192`, `heatwave.R:262-268` | The thin-climatology warning tests the maximum group size across cells. Two years of data give anomalies of exactly plus or minus 2 with no warning. `detrend = TRUE` skips the warning. (reproduced) | Test the minimum group size; apply in both branches. |
| 13 | `R/decompose.R:147-148` | A seasonal term is fitted with as few as one residual degree of freedom; per-cell slopes ranged from -0.14 to +0.80 degrees C per year on pure noise. No standard error is reported. (reproduced) | Require a residual-df margin; return a slope standard error. |
| 14 | `R/decompose.R:151,165` | A calendar month wholly `NA` for a cell makes `predict()` fail with "new levels". (reproduced) | Drop unused levels before prediction. |
| 15 | `R/temporal.R:230-240` | `integrate_covariate` adds `NA` as zero. A cell missing January to May returned 10 against 60. `methods.md` states the opposite. (reproduced) | Return `NA` below a minimum count, or add a count column. |
| 16 | `R/stratification.R:27-32, 74-77, 224` | The surface-to-bottom route is advertised, but `depths` takes two scalars, so N^2 and Eady shear are wrong wherever water depth differs. | Accept a per-cell depth column. |
| 17 | `R/watermass.R:76-77, 93-103` | Endmember names are ignored (swapped order gave 1.00 against 0.53); clamping and the residual are silent; example endmembers are uncited. (reproduced) | Validate names; warn on clamped share. |
| 18 | `R/ftle.R:126`, `fsle.R:127`, `checks.R:51` | Lagrangian warnings fire only at 90% `NA`. | Report the usable fraction. |
| 19 | `R/eke.R:105-108, 158-170` | EKE is exactly zero for one step or fewer than two years under `"climatology"`; rolling EKE is biased toward zero at record ends. (reproduced) | Warn on degenerate reference; document edge behaviour. |
| 20 | `R/fsle.R:231, 120` | Separation time is quantised to the step, biasing the exponent low (about 3% at 9 days, about 25% at 1 day for a 6 h step). | Interpolate the crossing time. |
| 21 | `R/residence.R:127-128, 178` | Censored values are recorded as `max_days` with no flag; the warning fires only above 10%. | Add a censoring indicator column. |
| 22 | `R/fronts.R:87-101` | A quantile threshold always produces fronts; a uniform-gradient field scored about 9% of cells as front. (reproduced) | Add a minimum absolute gradient. |
| 23 | `R/eddies.R:248, 34` | Polarity is wrong in the Southern Hemisphere. (reproduced) | Multiply by `sign(lat)`. |
| 24 | `R/transport.R:272-283, 145-146`; `methods.md:1045` | Section documentation is stale: the superseded score pair, normals of 259 degrees for both sections (Channel is 290.3), inconsistent depth profiles. | Update from current code. |
| 25 | `lcr-pipeline/` | Defaults use the surface level that the README says fails silently; `fetch.py:34-42` reuses cached files after a box change; the output column is `index`, not `LCR_extended`, and nothing in the file marks it as unvalidated; `run_record.py` records no parameters or versions. | Change defaults; key cache files by box; write a header and metadata. |
| 26 | `docs/lcr-extension-experiment.md` | The detrended correlation (r = -0.334, confidence interval excluding zero) is in the results table but not discussed; the headline is r = +0.18. | Discuss it beside the headline. |
| 27 | `.github/workflows/check-citations.yaml:18-24`; vignette | Path triggers omit `vignettes/**`, `man/**`, `NEWS.md`, `lcr-pipeline/**`; the reusable workflow is pinned to `@main`; the documented `inst/scripts/check_citations.R` does not exist; Ramp 1985 and Townsend 2015 have no registry rows. | Extend triggers; pin; add the script or correct the docs. |

### 5.3 Low severity

| # | Location | Defect |
|---|---|---|
| 28 | `R/deformation.R:99-114` | Spherical curvature terms omitted (about 7e-8 per second at 43 N, 0.5 m/s). |
| 29 | `R/gradient.R:74-75` | A diagonal `NA` neighbour removes a cell, widening the coastal `NA` band beyond the documented outer ring. |
| 30 | `R/gradient.R:104` | 111320 m per degree is the equatorial value, not the mean-radius value as documented; error about 0.1%. |
| 31 | `R/heatwave.R:138, 145` | Cold-spell `intensity` is negative and `cumulative` positive on the same event. |
| 32 | `R/temporal.R:234` | `per = "month"` divides by days / 30.4375; monthly steps give a plus or minus 5% error. |
| 33 | `R/shore.R:25-26` | Documented as ellipsoidal; `s2` computes on a sphere. Distance is unsigned. |
| 34 | `R/ftle.R:296-297` | Eigenvalues below 1 are clamped to zero; documented in `methods.md` but not in roxygen. |
| 35 | `R/stratification.R:171-190` | No CRS check on the latitude used by `eady_growth_rate`; a projected CRS gives plausible values. |
| 36 | `inst/CITATION` | `year` is `Sys.Date()`, so it changes with each call. |
| 37 | EML scale | `measurementScale` is `ratio` for all numeric columns, including Celsius, flags and categories. |

## 6. Limitations

- Literature claims were not checked against the papers. Items such as the
  Silver et al. (2023) figures, the Townsend et al. (2015) water-mass
  characterisation and the Grodsky et al. (2025) box are unverified.
- Section transports were not evaluated on real GLORYS data; the extracts that
  `docs/section-*.R` require were not present.
- The Python pipeline was not executed; its quoted statistics were not
  re-derived.
- The existing test suite was not run.
- Findings marked "reproduced" rest on one reviewer's snippet unless also
  marked "confirmed".

## 7. Verification performed

| Finding | Check |
|---|---|
| 1 | Code read; summation at `transport.R:101` treats dropped points as zero |
| 2 | Code read at `ftle.R:216-217` |
| 3 | Code read; run on a 60-day daily series, value 31 on days 1 and 10 |
| 4 | Code read at `region.R:65-71` |
| 5, 6, 8 | Code read in `eml-registry.R` and `eml.R` |

## 8. Planned work

Items 1 to 4 (silent wrong values) were corrected after this audit, each with a
regression test that fails on the previous code. See `NEWS.md`. Remaining order:
items 5 to 8 (archival metadata), then the warning-threshold changes in
section 4. Each fix should add a test that reproduces the defect first.

The item 2 correction has one consequence for users: the last steps of a
record, within `integration_days` of its end, now return `NA` for `ftle()`.

## 9. References

Belkin and O'Reilly (2009); d'Ovidio et al. (2004); Eady (1949); Haller
(2015); Hobday et al. (2016, 2018); Isern-Fontanet et al. (2003); Lindzen and
Farrell (1980); Okubo (1970); Ramp et al. (1985); Townsend et al. (2015);
UNESCO (1983); Weiss (1991). Full entries are in
`vignettes/derivoce.Rmd`.
