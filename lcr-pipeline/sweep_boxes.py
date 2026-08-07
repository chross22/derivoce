"""Ask whether any arrival-box placement makes the index track the published series.

STATUS.md has named the boxes as the most likely thing to be wrong from the
start: they stand in for hydrographic sections whose coordinates the paper does
not publish, and the Labrador box in particular sits so close to the release
line that it measures "drifted a bit south" rather than a crossing. Where the
boxes sit does not change the trajectories, only the counting, so this costs a
re-tally rather than a re-advection.

Read the whole distribution, not the best cell. Searching a grid of boxes
against a target series is fitting: with enough candidates one will correlate
by chance, and the honest question is whether the good placements are a
coherent region or a scattering of lucky cells. A placement chosen here has to
be confirmed on years this sweep never saw.

    python sweep_boxes.py --depth 47.374
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from validate import PUBLISHED, deseasonalize

HERE = Path(__file__).parent
OUT = HERE / "output"

# Candidate boxes. Labrador is moved south, away from the release line; Scotian
# is moved along the shelf. Both keep the paper's shape and only slide.
LABRADOR_BANDS = [(50, 52), (49, 51), (48, 50), (47, 49), (46, 48), (45, 47)]
SCOTIAN_BANDS = [(-66, -59), (-64, -57), (-62, -55), (-60, -53), (-58, -51)]
LABRADOR_LON = (-56.0, -52.0)
SCOTIAN_LAT = (42.5, 45.5)


def load_blocks(depth: float) -> tuple:
    """Every block's trajectories, with each particle's release date."""
    zarrs = sorted(OUT.glob(f"traj_*_{depth:g}m.zarr"))
    if not zarrs:
        raise SystemExit(f"no trajectories at {depth:g} m in {OUT}")

    lons, lats, dates = [], [], []
    for z in zarrs:
        csv = OUT / z.name.replace("traj_", "lcr_releases_").replace(".zarr", ".csv")
        if not csv.exists():
            print(f"  skipping {z.name}: no matching release list")
            continue
        released = pd.read_csv(csv)["DATE"].to_numpy()
        t = xr.open_zarr(z)
        # zarr round-trips the write-once release id as float, since the
        # variable carries a NaN fill value.
        which = np.asarray(t["release"].values).reshape(-1).astype(np.int64)
        lons.append(t["lon"].values)
        lats.append(t["lat"].values)
        dates.append(released[which])
        print(f"  {z.name}: {len(which)} particles, {len(released)} releases")

    width = max(a.shape[1] for a in lons)
    def pad(a):
        out = np.full((a.shape[0], width), np.nan)
        out[:, :a.shape[1]] = a
        return out

    return (np.vstack([pad(a) for a in lons]),
            np.vstack([pad(a) for a in lats]),
            pd.to_datetime(np.concatenate(dates)))


def inside(lon, lat, box) -> np.ndarray:
    return ((lon >= box["lon"][0]) & (lon <= box["lon"][1])
            & (lat >= box["lat"][0]) & (lat <= box["lat"][1])).any(axis=1)


def monthly(dates: pd.DatetimeIndex, lab: np.ndarray, sco: np.ndarray) -> pd.DataFrame:
    df = pd.DataFrame({"YEAR": dates.year, "MONTH": dates.month,
                       "lab": lab, "sco": sco})
    g = df.groupby(["YEAR", "MONTH"], as_index=False).agg(
        released=("lab", "size"), labrador=("lab", "sum"), scotian=("sco", "sum"))
    g["INDEX"] = (g["labrador"] - g["scotian"]) / g["released"]
    return g


def score(g: pd.DataFrame, pub: pd.DataFrame) -> tuple:
    m = pub.merge(g, on=["YEAR", "MONTH"])
    if len(m) < 24 or m["YEAR"].nunique() < 5:
        return np.nan, np.nan
    a = deseasonalize(m, "LCR")
    b = deseasonalize(m, "INDEX")
    if not np.isfinite(a).all() or np.std(b) == 0:
        return np.nan, np.nan
    return float(np.corrcoef(a, b)[0, 1]), len(m)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--depth", type=float, default=47.374)
    args = p.parse_args()

    print("loading trajectories")
    lon, lat, dates = load_blocks(args.depth)
    print(f"  {lon.shape[0]} particles total, "
          f"{dates.min():%Y-%m} to {dates.max():%Y-%m}\n")

    pub = pd.read_csv(PUBLISHED)
    pub.columns = [c.upper() for c in pub.columns]

    # Each candidate box is tested against the trajectories once, then reused
    # across every pairing.
    lab_hits = {band: inside(lon, lat, dict(lon=LABRADOR_LON, lat=band))
                for band in LABRADOR_BANDS}
    sco_hits = {band: inside(lon, lat, dict(lon=band, lat=SCOTIAN_LAT))
                for band in SCOTIAN_BANDS}

    print("deseasonalized r against the published series")
    print("rows: Labrador box latitude   cols: Scotian box longitude\n")
    header = "  LAB \\ SCO " + "".join(f"{f'{b[0]},{b[1]}':>12}" for b in SCOTIAN_BANDS)
    print(header)

    results = []
    for lband in LABRADOR_BANDS:
        cells = []
        for sband in SCOTIAN_BANDS:
            g = monthly(dates, lab_hits[lband], sco_hits[sband])
            r, n = score(g, pub)
            results.append((r, lband, sband, n))
            cells.append(f"{r:+12.3f}" if np.isfinite(r) else f"{'--':>12}")
        star = " *" if lband == (50, 52) else "  "
        print(f"  {lband[0]}-{lband[1]}N{star}" + "".join(cells))
    print("\n  * the box as currently configured")

    good = [x for x in results if np.isfinite(x[0])]
    if not good:
        raise SystemExit("\nno placement produced a usable series")

    rs = np.array([x[0] for x in good])
    best = max(good, key=lambda x: x[0])
    print(f"\n  {len(good)} placements tested")
    print(f"  r ranges {rs.min():+.3f} to {rs.max():+.3f}, median {np.median(rs):+.3f}")
    print(f"  best: Labrador {best[1][0]}-{best[1][1]}N, "
          f"Scotian {best[2][0]},{best[2][1]}  r {best[0]:+.3f} on {best[3]} months")
    print(f"  placements above 0.4: {(rs > 0.4).sum()} of {len(rs)}")

    if rs.max() < 0.4:
        print("\n  No placement comes close. The disagreement is not the boxes.")
    else:
        print("\n  Chosen by searching against the target, so this r is optimistic.")
        print("  Confirm it on years this sweep never saw before believing it.")


if __name__ == "__main__":
    main()
