# Where this was left off

**Written, never run.** Every script parses and its `--help` works. No trajectory
has been computed and nothing has been validated. Treat all of it as untested.

## What exists

| File | State |
|---|---|
| `README.md` | Method, and where it departs from the paper |
| `fetch.py` | Written, not run. Downloads daily GLORYS, one file per year, resumable |
| `track.py` | Written, not run. Seeds, advects 3 years, counts arrivals |
| `validate.py` | Written, not run. Correlates against the published series |
| `data/lcr_published.csv` | **Ready.** 260 months, 1993-2014, exported from datamatch |

`parcels 3.1.4` is installed in `~/miniconda3` and imports cleanly.

## The next three commands

```bash
cd lcr-pipeline

# 1. Two years is enough to find out whether any of this works. ~250 MB.
python fetch.py --start 1995-01-01 --end 1999-12-31

# 2. One year of monthly releases. Trajectories run 3 years past the last
#    release, which is why step 1 fetches to 1999.
python track.py --start 1995-01-01 --end 1995-12-31 --particles 100

# 3. Twelve months is thin, but a correlation near zero is still informative.
python validate.py
```

Start small. A release of 100 particles tracked three years at a six-hour step is
not instant, and finding out that the arrival regions are wrong is much cheaper
at twelve releases than at two hundred.

## What is most likely to be wrong

**The arrival regions.** `LABRADOR` and `SCOTIAN` in `track.py` are boxes
standing in for hydrographic sections whose coordinates the paper does not
publish. This is the first thing to change if `validate.py` reports a poor
correlation, and it is worth plotting a few trajectories before assuming the
physics is at fault.

**`FieldSet.from_data` loading everything into memory.** Fine for a few years,
not for the full 1993-present record at ~4 GB. If it becomes a problem, switch to
`FieldSet.from_netcdf` with a file list, which streams.

**Whether daily is enough.** The monthly attempt failed because averaging removed
the Labrador Current jet and the Grand Banks bifurcation
(`../docs/lcr-extension-experiment.md`). Daily should retain both, since that is
what Jutras et al. used, but it is an assumption until step 3 says otherwise.

## What would count as success

`validate.py` reporting a correlation above about 0.7 on the 1993-2014 overlap.
The absolute values will not match, and are not meant to: this counts particles
entering boxes, theirs counts crossings of sections. Whether the two move
together is the question.

If it does track, the extension past 2014 is worth computing and can be joined as
a covariate on `YEAR`/`MONTH` like any other monthly index. If it does not, say
so in `../docs/lcr-extension-experiment.md` alongside the monthly result, and the
honest conclusion is that the index cannot be reproduced outside the original
framework.

## What this is not

Not Jutras et al.'s series. Different integrator, different release cadence,
different arrival criterion. If anything from it is published, cite them for the
method and describe this as an independent recomputation.
