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

## The domain: tested, and it did not matter

`BOX` was extended east from 45W to 30W and south from 40N to 36N, and the 1995
pilot rerun on both. **The index barely moved.** Monthly Scotian counts differ
by 0-3 particles, the largest change in the index is 0.006, and the two series
correlate at 1.000.

| | median life | completed | LAB | SCO |
|---|---|---|---|---|
| 45W edge | 600 d | 39% | 55.8% | 25.3% |
| 30W edge | 900 d | 44% | 55.8% | 25.5% |

This refuted a plausible-looking argument, which is worth recording so it is not
re-derived. Losses looked like a bias: Labrador arrivals complete within 165
days, Scotian arrivals need 410 (median), and particles reaching neither box
left at a median of 130 days -- so the Labrador term appeared to be counted in
full while the Scotian term was truncated. Empirically it is not. Given 15 more
degrees of ocean, those particles still reach neither box; they are being
exported, not cut off.

Two things did improve, and are why the larger box is kept: retention rose from
600 to 900 days, and the Scotian-band distribution now **closes** inside the
domain (counts decay to 5.8% at the eastern edge) instead of still climbing at
the boundary. A distribution that peaks at the edge is measuring the edge.

The cost is real: 186 MB per year against 95 MB, so about 4.6 GB for 1993-2017,
and ~1.5 GB per worker. If that becomes the binding constraint, the smaller box
is defensible on this evidence -- it gives the same answer.

## What is still open

**The Labrador box sits next to the release line.** Median time-to-arrival is
30 days, and sweeping the box south gives a smooth distance-decay curve
(75% → 4%) rather than a plateau. It is measuring "drifted a bit south", not a
section crossing, and it carries most of the index's magnitude. Moving it well
downstream is the next methodological change worth testing, and `tally_releases`
takes both regions as arguments so a sweep costs a re-tally, not a re-advection.

**Nothing is validated yet.** Only 1995 has been computed, and `validate.py`
withholds both the deseasonalized correlation and any verdict on one year, on
purpose.

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

`run_record.py` does this and has not yet been run in anger.

Releases in year Y need velocities for Y..Y+3 and nothing else, so each year is
an independent job. It runs them as separate processes, which keeps every
fieldset to one four-year window rather than the whole record, and merges the
per-release output at the end.

```bash
python fetch.py --start 1993-01-01 --end 2017-12-31 --depth 47.374 --jobs 4
python run_record.py --start 1993 --end 2014 --depth 47.374 --jobs 4
python validate.py --depth 47.374
```

`--jobs` is limited by memory, not cores: each worker holds four years of daily
velocities, about 1.5 GB at the current box. Four workers is roughly 6 GB on a
34 GB machine. Separate processes rather than a thread pool because parcels
compiles a kernel per ParticleSet and writes C to a scratch directory.

Budget: ~4.6 GB of downloads for 1993-2017, and 22 blocks of ~5300 particles.

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
