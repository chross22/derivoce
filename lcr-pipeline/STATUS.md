# Where this was left off

The pipeline **runs**, is **tested**, and has been **exercised on 1995**. It has
not yet been validated against the published series: that needs several years of
overlap, and only one year has been computed.

The 1995 pilots settled two questions the previous pass could not: which depth
to advect on, and whether the domain is large enough. The answers are 47.374 m
and no.

## What exists

| File | State |
|---|---|
| `README.md` | Method, and where it departs from the paper |
| `fetch.py` | Runs. Daily GLORYS, one file per year per depth, parallel across years, resumable |
| `track.py` | Runs. Seeds every release into one ParticleSet and advects them in a single pass |
| `diagnose.py` | Runs. Survival, exit boundary, occupancy, and box sweeps for one configuration |
| `compare.py` | Runs. Several configurations side by side |
| `validate.py` | Runs. Deseasonalized correlation against the published series, with autocorrelation-adjusted significance |
| `tests/` | 44 tests, all passing (`pytest tests/`) |
| `data/lcr_published.csv` | 260 months, 1993-2014, exported from datamatch |

`parcels 3.1.4`, `copernicusmarine 2.2.3`, `pytest 9.1.1` in `~/miniconda3`.

Velocities on disk: 1995-1998 at 0.494, 47.374, 55.764 and 92.326 m (~1.5 GB).
All gitignored.

## What the pilots showed

One year of weekly releases, 100 particles each, 5300 particles per run.

| depth | median life | completed 3 yr | LAB arrivals | SCO arrivals | SCO box masked |
|---|---|---|---|---|---|
| 0.494 m | 155 d | 24% | 52.1% | **1.1%** | 19.4% |
| **47.374 m** | 600 d | 39% | 55.8% | **25.3%** | 29.6% |
| 55.764 m | 555 d | 39% | 55.9% | 19.9% | 32.7% |
| 92.326 m | 1115 d | 54% | 53.5% | 5.6% | 46.3% |

**Use 47.374 m.** At the surface, Ekman drift sweeps particles offshore and
almost nothing reaches the Scotian Shelf, so the index is the Labrador count
alone. Below about 65 m the shelf is *land* -- GLORYS masks each level against
the bathymetry, and at 92 m 57% of the Grand Banks and 46% of the Scotian box
are dry, so particles stall against the mask instead of crossing. 47 m sits
below the Ekman layer and above the shelf. `track.py` now prints the masked
fraction of both boxes on every run; do not ignore it.

## What is still wrong

**The domain is too small on the east.** `BOX` in `fetch.py` stops at 45W. At
47 m, 3246 of 5300 particles leave early and 2165 of those cross the eastern
edge -- **72% of them at or below 48N**, peaking at 44-46N. That is the Grand
Banks slope water, not the North Atlantic Current: particles that could still
have turned southwest onto the shelf.

This biases the index rather than merely thinning it:

- Labrador arrivals are complete -- 99% of them happen within **165 days**.
- Scotian arrivals need **410 days** (median), 620 at the 90th percentile.
- Particles that reached neither box and left early lasted a median of **130 days**.

So the Labrador term is counted in full and the Scotian term is truncated before
most of its arrivals could occur, by an amount that varies month to month with
the flow. Fix the domain before reading anything into a correlation.

Suggested: extend `BOX` east to about 30W and south to about 35N, and rerun the
pilot before committing to the full record.

**Also unresolved:** the Labrador box sits next to the release line. Median
time-to-arrival is 30 days, and sweeping the box south gives a smooth
distance-decay curve (75% → 4%). It is measuring "drifted a bit south", not a
section crossing. Consider moving it well downstream once the domain is fixed.

## Running it

```bash
python fetch.py --start 1995-01-01 --end 1998-12-31 --depth 47.374 --jobs 4
python track.py --start 1995-01-01 --end 1995-12-31 --particles 100 --depth 47.374
python diagnose.py --depth 47.374
python validate.py --depth 47.374
```

`track.py` needs fields to three years past the last release, so a run of
releases in year Y wants Y through Y+3 fetched.

## Scaling to the full record

Not yet attempted. Two known obstacles:

1. `FieldSet.from_data` holds every field in memory. Four years at the current
   box is ~760 MB; the full 1993-2017 record at an enlarged box would be ~10 GB.
   Either switch to `FieldSet.from_netcdf`, which streams, or run in **blocks**:
   releases in year Y need only Y..Y+3, so each year is an independent job. That
   parallelises across cores and keeps each fieldset small. Blocks are the
   cheaper change and the machine has 10 cores.
2. Download is ~185 MB per year at the enlarged box, so ~4.6 GB for 1993-2017.

## What would count as success

`validate.py` reporting a **deseasonalized** correlation above about 0.7 over
the 1993-2014 overlap. Not the raw correlation: both series carry the annual
cycle of the Labrador Current, and two series sharing a seasonal cycle correlate
highly whether or not they measure the same thing. `validate.py` withholds the
deseasonalized number below five years of overlap and withholds any verdict
below 24 months, on purpose.

Absolute values will not match and are not meant to: this counts particles
entering boxes, theirs counts crossings of sections.

If it does track, the extension past 2014 can be joined as a covariate on
`YEAR`/`MONTH` like any other monthly index. If it does not, say so in
`../docs/lcr-extension-experiment.md` alongside the monthly result, and the
honest conclusion is that the index cannot be reproduced outside the original
framework.

## What this is not

Not Jutras et al.'s series. Different integrator, different release cadence,
different arrival criterion. If anything from it is published, cite them for the
method and describe this as an independent recomputation.
