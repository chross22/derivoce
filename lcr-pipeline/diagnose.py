"""Ask where the particles actually went, before asking whether the index is right.

A correlation against the published series says a configuration is wrong but not
why. This says why: how many particles survived, which boundary they left
through, where the cloud sat, and how the arrival counts respond to moving the
boxes. Run it on every configuration before trusting a correlation from it.

    python diagnose.py --depth 0.494
    python diagnose.py --depth 92.326 --start 199501 --end 199512
"""
import argparse
from pathlib import Path

import numpy as np
import xarray as xr

from fetch import BOX, SURFACE
from track import LABRADOR, SCOTIAN, RELEASE, TRACK_DAYS

OUT = Path(__file__).parent / "output"


def load(depth: float, start: str, end: str) -> xr.Dataset:
    path = OUT / f"traj_{start}_{end}_{depth:g}m.zarr"
    if not path.exists():
        raise SystemExit(f"no trajectories at {path}. Run track.py --depth {depth:g}")
    print(f"{path.name}\n")
    return xr.open_zarr(path)


def survival(lon: np.ndarray, n: int, n_obs: int) -> np.ndarray:
    """How much of the intended three years the ensemble actually got."""
    valid = np.isfinite(lon)
    lifespan = valid.sum(axis=1) * 5
    full = (lifespan >= TRACK_DAYS).sum()

    print("SURVIVAL")
    print(f"  particles            {n}")
    print(f"  median lifespan      {np.median(lifespan):.0f} of {TRACK_DAYS} days")
    print(f"  completed full track {full} ({100 * full / n:.0f}%)")
    if full < 0.5 * n:
        print("  ^ most of the ensemble is truncated; counts below are not "
              "comparable across releases")
    return valid


def exits(lon: np.ndarray, lat: np.ndarray, valid: np.ndarray, n: int) -> None:
    """Which boundary the lost particles crossed.

    Losses concentrated on one edge mean the domain is too small there, not
    that the particles went anywhere meaningful.
    """
    lifespan = valid.sum(axis=1) * 5
    died = lifespan < TRACK_DAYS
    if not died.any():
        print("\nEXITS\n  none: every particle completed the track")
        return

    n_obs = lon.shape[1]
    last = np.where(valid[died], np.arange(n_obs), -1).max(axis=1)
    idx = np.arange(died.sum())
    final_lon = lon[died][idx, last]
    final_lat = lat[died][idx, last]

    edge = 0.5
    on = {
        "west  ": final_lon <= BOX["lon"][0] + edge,
        "east  ": final_lon >= BOX["lon"][1] - edge,
        "south ": final_lat <= BOX["lat"][0] + edge,
        "north ": final_lat >= BOX["lat"][1] - edge,
    }
    print(f"\nEXITS  ({died.sum()} of {n} particles lost early)")
    for name, mask in on.items():
        k = int(mask.sum())
        print(f"  {name} {k:6d}  ({100 * k / died.sum():5.1f}% of losses)"
              f"  {'#' * int(40 * k / died.sum())}")
    interior = int((~np.logical_or.reduce(list(on.values()))).sum())
    print(f"  interior {interior:4d}  ({100 * interior / died.sum():5.1f}%)"
          "   grounded or hit a land mask")


def occupancy(lon: np.ndarray, lat: np.ndarray) -> None:
    ok = np.isfinite(lon) & np.isfinite(lat)
    H, ye, xe = np.histogram2d(
        lat[ok], lon[ok],
        bins=[np.arange(BOX["lat"][0], BOX["lat"][1] + 2, 2),
              np.arange(BOX["lon"][0], BOX["lon"][1] + 2, 2)])
    H = 100 * H / H.sum()
    print("\nOCCUPANCY  (% of all valid positions, 2-degree cells)")
    print("       " + "".join(f"{x:6.0f}" for x in xe[:-1]))
    for i in range(len(ye) - 1)[::-1]:
        row = "".join(f"{v:6.1f}" if v >= 0.05 else "     ." for v in H[i])
        print(f"  {ye[i]:3.0f}N {row}")


def hits(lon, lat, lo0, lo1, la0, la1) -> int:
    inside = (lon >= lo0) & (lon <= lo1) & (lat >= la0) & (lat <= la1)
    return int(inside.any(axis=1).sum())


def sweep(lon: np.ndarray, lat: np.ndarray, n: int) -> None:
    """Move the boxes and watch the counts.

    Trajectories do not depend on where the boxes are, so this is free. A count
    that keeps climbing towards a domain edge is measuring the edge; a count
    that decays smoothly with distance from the release line is measuring
    distance, not a crossing.
    """
    print("\nSWEEP: Scotian band (lat 42.5-45.5), west to east")
    for lo0 in range(int(BOX["lon"][0]), int(BOX["lon"][1]) - 1, 2):
        h = hits(lon, lat, lo0, lo0 + 2, 42.5, 45.5)
        flag = "  <- SCOTIAN" if lo0 >= SCOTIAN["lon"][0] and lo0 < SCOTIAN["lon"][1] else ""
        print(f"  {lo0:4d} to {lo0+2:4d}  {h:6d}  {100*h/n:5.1f}%  "
              f"{'#' * int(40 * h / n)}{flag}")

    print("\nSWEEP: Labrador box moved south from the release line "
          f"(lat {RELEASE[0][1]:.0f}-{RELEASE[1][1]:.1f})")
    for la0 in np.arange(52, 44, -1.0):
        h = hits(lon, lat, LABRADOR["lon"][0], LABRADOR["lon"][1], la0, la0 + 2)
        flag = "  <- LABRADOR" if la0 == LABRADOR["lat"][0] else ""
        print(f"  lat {la0:4.0f} to {la0+2:4.0f}  {h:6d}  {100*h/n:5.1f}%  "
              f"{'#' * int(40 * h / n)}{flag}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--depth", type=float, default=SURFACE)
    p.add_argument("--start", default="199501", help="YYYYMM of first release")
    p.add_argument("--end", default="199512", help="YYYYMM of last release")
    args = p.parse_args()

    traj = load(args.depth, args.start, args.end)
    lon, lat = traj["lon"].values, traj["lat"].values
    n, n_obs = lon.shape

    valid = survival(lon, n, n_obs)
    exits(lon, lat, valid, n)
    occupancy(lon, lat)

    print(f"\nARRIVALS with the boxes as configured")
    print(f"  LABRADOR {hits(lon, lat, *LABRADOR['lon'], *LABRADOR['lat']):6d}"
          f"  ({100*hits(lon, lat, *LABRADOR['lon'], *LABRADOR['lat'])/n:5.1f}%)")
    print(f"  SCOTIAN  {hits(lon, lat, *SCOTIAN['lon'], *SCOTIAN['lat']):6d}"
          f"  ({100*hits(lon, lat, *SCOTIAN['lon'], *SCOTIAN['lat'])/n:5.1f}%)")

    sweep(lon, lat, n)


if __name__ == "__main__":
    main()
