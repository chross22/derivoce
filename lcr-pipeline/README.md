# Extending the Labrador Current retroflection index

`LCR` in datamatch is Jutras et al. (2023) and stops at 2014. This recomputes it
from Copernicus daily fields with [OceanParcels](https://oceanparcels.org/), so
it can run to the present.

**Separate from the R package on purpose.** derivoce needs only `sf` and `terra`.
Making it depend on a Python scientific stack would break `R CMD check` anywhere
that stack is missing. The boundary between the two is a CSV: this writes
`YEAR,MONTH,LCR_extended`, which joins on date like any other monthly index.

Why Python at all: an earlier attempt reused derivoce's own RK4 advection on
monthly fields and failed, for reasons recorded in
[`../docs/lcr-extension-experiment.md`](../docs/lcr-extension-experiment.md).
Monthly averaging removes the narrow Labrador Current jet and the Grand Banks
bifurcation the index depends on. Daily fields are needed, and at ~1 billion
interpolations that is not an R job.

## Status

Pipeline written, **not yet validated**. Nothing here should be trusted until
`validate.py` shows it tracks the published series over 1993-2014.

## Setup

```bash
pip install parcels          # brings xarray, netCDF4, zarr, scipy
copernicusmarine login       # once
```

## Running

```bash
# 1. Daily velocities, fetched in parallel across years.
python fetch.py --start 1995-01-01 --end 1998-12-31 --depth 47.374 --jobs 4

# 2. Track particles and count arrivals
python track.py --start 1995-01-01 --end 1995-12-31 --particles 100 --depth 47.374

# 3. Where did the particles actually go?
python diagnose.py --depth 47.374

# 4. Compare against the published series over the overlap
python validate.py --depth 47.374
```

`fetch.py` writes one NetCDF per year per depth to `data/`. Trajectories run
three years past the last release, so releases in year Y need Y through Y+3.

Step 3 is not optional. A correlation says a configuration is wrong but not why;
`diagnose.py` reports particle survival, which boundary the losses crossed, and
how the arrival counts respond to moving the boxes. Run `compare.py` to put
several configurations side by side.

```bash
pytest tests/
```

## The method, and where it departs from the paper

Jutras et al. seed particles weekly across (53N, 56.7W)-(54.3N, 52.0W), track
each for three years through GLORYS12V1, and take the difference between counts
reaching the Labrador and Scotian shelves. High retroflection means fewer
particles continue southwest.

This follows that, and differs in ways that matter:

- **Release cadence and count** are arguments here, not fixed at their 966 per
  week. Fewer particles is faster and noisier.
- **Arrival is counted by region**, not by crossing a named hydrographic section.
  Their sections are not published as coordinates, so boxes stand in. This is the
  largest source of disagreement and the first thing to revisit if validation is
  poor.
- **Advection on one level, at 47.374 m.** Depth turned out to matter more than
  anything else. At the surface, Ekman drift carries particles offshore and
  almost none reach the Scotian Shelf. Below about 65 m the shelf is land:
  GLORYS masks each level against the bathymetry, so at 92 m most of the Grand
  Banks is dry and particles stall against the mask. 47 m is below the Ekman
  layer and above the shelf. See `STATUS.md` for the comparison.

So this is a reimplementation, not their series. It has to be validated on the
overlap before the extension means anything, and if published, cite them for the
method and describe this as an independent recomputation.

## Reference

Jutras M, Dufour CO, Mucci A, Talbot LC (2023). Large-scale control of the
retroflection of the Labrador Current. *Nature Communications* **14**, 2623.
[doi:10.1038/s41467-023-38321-y](https://doi.org/10.1038/s41467-023-38321-y)
