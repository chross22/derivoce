# Where this was left off

**The pipeline works and the answer is negative.** It runs, it is tested, and it
produces a stable index over 1993-1997 that does **not** track the published
series: deseasonalized r = +0.18 (p = 0.36, 95% CI [-0.13, +0.28]) on 56 months.

The full reasoning is in
[`../docs/lcr-extension-experiment.md`](../docs/lcr-extension-experiment.md),
alongside the earlier monthly-fields failure. Read that first. This file is the
operational detail.

Do not attach anything from this as a covariate.

## Why it is not a bug

Three explanations were tested and eliminated, which is why the conclusion is
about the method rather than the implementation.

- **Not the boxes.** 30 placements swept over the same 26,500 particles.
  Deseasonalized r ranges -0.162 to +0.243, median -0.050, nothing above 0.4.
  The configured box is already near the best. `python sweep_boxes.py`
- **Not the sample size.** Split-half reliability is 0.827, so the ceiling on
  any correlation is 0.91. The index is highly reproducible; more particles
  would not help. `python reliability.py`
- **Not the domain.** Extending the eastern edge from 45W to 30W raised
  retention from 600 to 900 days and moved the index by at most 0.006.

What remains is the one departure that could not be tested: this counts
particles **entering a region**, the paper counts particles **crossing a
hydrographic section**, and the section coordinates are not published.

## What exists

| File | State |
|---|---|
| `README.md` | Method, and where it departs from the paper |
| `fetch.py` | Daily GLORYS, one file per year per depth, parallel across years, atomic and resumable |
| `track.py` | Seeds every release into one ParticleSet and advects them in a single pass |
| `run_record.py` | Runs the record as parallel one-year blocks and merges |
| `diagnose.py` | Survival, exit boundary, occupancy, box sweeps for one configuration |
| `compare.py` | Several configurations side by side |
| `sweep_boxes.py` | Arrival-box placements scored against the published series |
| `reliability.py` | Split-half reliability: signal against counting noise |
| `validate.py` | Deseasonalized correlation, autocorrelation-adjusted significance |
| `tests/` | 50 tests (`pytest tests/`) |
| `data/lcr_published.csv` | 260 months, 1993-2014, exported from datamatch |

`parcels 3.1.4`, `copernicusmarine 2.2.3`, `pytest 9.1.1` in `~/miniconda3`.

On disk: 1993-2014 at 47.374 m (~4 GB), plus 1995-1998 at 0.494, 55.764 and
92.326 m under `data/small/` at the old smaller domain. All gitignored.
Trajectories for 1993-1997 are in `output/`.

## Use 47.374 m

Depth is the one setting that decides whether this works at all, and the window
is narrow. GLORYS masks each level against the bathymetry, so a level deeper
than the shelf turns the shelf into land.

| level | median lifespan | completed | LAB | SCO | SCO box masked |
|---|---|---|---|---|---|
| 0.494 m | 155 d | 24% | 52.1% | **1.1%** | 19.4% |
| **47.374 m** | 600 d | 39% | 55.8% | **25.3%** | 29.6% |
| 55.764 m | 555 d | 39% | 55.9% | 19.9% | 32.7% |
| 92.326 m | 1115 d | 54% | 53.5% | 5.6% | 46.3% |

At the surface Ekman drift carries particles offshore; below about 65 m the
shelf is dry. A run at the wrong level fails **silently**, producing near-zero
counts indistinguishable from a real absence of transport, so `track.py` prints
the masked fraction of both boxes on every run. Do not ignore it.

Depths must be actual GLORYS levels. `fetch.py` asks for a narrow band around
the value and a band falling between levels returns nothing.

## Running it

```bash
python fetch.py --start 1993-01-01 --end 2017-12-31 --depth 47.374 --jobs 4
python run_record.py --start 1993 --end 2014 --depth 47.374 --jobs 4
python validate.py --depth 47.374
```

Releases in year Y need velocities for Y..Y+3, so releases through 2014 need
fields to 2017. `--jobs` is limited by memory, not cores: each worker holds four
years of daily velocities, about 1.5 GB, so four workers is roughly 6 GB.

For a single year, `track.py` takes the same arguments and writes output tagged
with the release span. Always run `diagnose.py` before reading a correlation --
it reports survival, which boundary the losses crossed, and how the counts
respond to moving the boxes. A correlation says a configuration is wrong; only
the diagnostics say why.

## If this is picked up again

The one thing that would change the answer is the **section coordinates** from
Jutras et al. Everything else has been tested. With them, replace the box test
in `tally_releases` with a crossing test; it takes both regions as arguments
precisely so the counting rule can be changed without re-advecting.

Two smaller things, if the definition question is ever resolved:

- The Labrador box sits next to the release line. Median time-to-arrival is 30
  days and sweeping it south gives a smooth distance-decay curve, so it is
  largely measuring "drifted a bit south".
- Only 1993-1997 has been advected. 1993-2014 velocities are on disk, so the
  remaining blocks are a `run_record.py` call away.

## What this is not

Not Jutras et al.'s series, and now demonstrably not equivalent to it. If
anything from this is published, cite them for the method and describe it as an
independent recomputation that did not reproduce their index.
