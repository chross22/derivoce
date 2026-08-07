"""Compare configurations side by side.

diagnose.py explains one run in depth; this puts several next to each other, so
a choice between them can be made on retention and reachability rather than on
a correlation that may be measuring the domain edge.

    python compare.py --depths 0.494 47.374 55.764 92.326
"""
import argparse
from pathlib import Path

import numpy as np
import xarray as xr

from fetch import BOX, velocity_file
from track import LABRADOR, SCOTIAN, TRACK_DAYS, masked_fraction

OUT = Path(__file__).parent / "output"


def hits(lon, lat, region) -> int:
    inside = ((lon >= region["lon"][0]) & (lon <= region["lon"][1])
              & (lat >= region["lat"][0]) & (lat <= region["lat"][1]))
    return int(inside.any(axis=1).sum())


def summarise(depth: float, start: str, end: str, year: int) -> dict | None:
    path = OUT / f"traj_{start}_{end}_{depth:g}m.zarr"
    if not path.exists():
        return None
    traj = xr.open_zarr(path)
    lon, lat = traj["lon"].values, traj["lat"].values
    n, n_obs = lon.shape

    valid = np.isfinite(lon)
    lifespan = valid.sum(axis=1) * 5
    died = lifespan < TRACK_DAYS

    east = np.nan
    if died.any():
        last = np.where(valid[died], np.arange(n_obs), -1).max(axis=1)
        final_lon = lon[died][np.arange(died.sum()), last]
        east = 100 * (final_lon >= BOX["lon"][1] - 0.5).mean()

    vel = velocity_file(year, depth)
    masks = {}
    if vel.exists():
        ds = xr.open_dataset(vel)
        if "depth" in ds.dims:
            ds = ds.isel(depth=0)
        masks = {"lab": 100 * masked_fraction(ds, LABRADOR),
                 "sco": 100 * masked_fraction(ds, SCOTIAN)}

    return dict(
        depth=depth,
        median_life=float(np.median(lifespan)),
        completed=100 * (lifespan >= TRACK_DAYS).mean(),
        east_exits=east,
        lab=100 * hits(lon, lat, LABRADOR) / n,
        sco=100 * hits(lon, lat, SCOTIAN) / n,
        lab_mask=masks.get("lab", np.nan),
        sco_mask=masks.get("sco", np.nan),
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--depths", type=float, nargs="+", required=True)
    p.add_argument("--start", default="199501")
    p.add_argument("--end", default="199512")
    args = p.parse_args()

    rows = [r for d in args.depths
            if (r := summarise(d, args.start, args.end, int(args.start[:4])))]
    if not rows:
        raise SystemExit("no trajectory files found for those depths")

    print(f"{'depth':>9} | {'median':>7} {'done':>6} {'E-exit':>7} | "
          f"{'LAB':>6} {'SCO':>6} | {'LABmask':>8} {'SCOmask':>8}")
    print(f"{'(m)':>9} | {'life(d)':>7} {'(%)':>6} {'(%loss)':>7} | "
          f"{'(%)':>6} {'(%)':>6} | {'(%)':>8} {'(%)':>8}")
    print("-" * 76)
    for r in rows:
        print(f"{r['depth']:9.3f} | {r['median_life']:7.0f} {r['completed']:6.0f} "
              f"{r['east_exits']:7.0f} | {r['lab']:6.1f} {r['sco']:6.1f} | "
              f"{r['lab_mask']:8.1f} {r['sco_mask']:8.1f}")

    print(f"\nTRACK_DAYS is {TRACK_DAYS}; median life below that means the "
          "ensemble is\ntruncated and releases are not comparable to each other.")
    print("A high SCOmask means the Scotian box is land at that depth, so a low\n"
          "SCO arrival rate says nothing about transport.")


if __name__ == "__main__":
    main()
