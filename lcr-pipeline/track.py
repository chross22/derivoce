"""Track particles from the Labrador Shelf and count where they arrive.

Seeds a line across the Labrador Shelf at intervals, advects each release for
three years through daily GLORYS velocities, and records how many reach the
Labrador Shelf and how many reach the Scotian Shelf.

The index is the difference between those counts, scaled by the number released.
High retroflection means fewer particles carry on southwest, so the difference
grows.
"""
import argparse
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from parcels import AdvectionRK4, FieldSet, JITParticle, ParticleSet

DATA = Path(__file__).parent / "data"
OUT = Path(__file__).parent / "output"

# Jutras et al.'s release line, (lon, lat) at each end.
RELEASE = ((-56.7, 53.0), (-52.0, 54.3))

# Arrival regions. These stand in for the paper's hydrographic sections, whose
# coordinates are not published, and are the most likely thing to need changing
# if validation is poor.
LABRADOR = dict(lon=(-56.0, -52.0), lat=(50.0, 52.0))
SCOTIAN = dict(lon=(-66.0, -59.0), lat=(42.5, 45.5))

TRACK_DAYS = 3 * 365


def build_fieldset(years: list[int]) -> FieldSet:
    """One fieldset spanning every year needed, so trajectories cross year ends."""
    files = [DATA / f"glorys_{y}.nc" for y in years]
    missing = [f.name for f in files if not f.exists()]
    if missing:
        raise SystemExit(f"missing velocity files: {missing}. Run fetch.py first.")

    ds = xr.open_mfdataset([str(f) for f in files], combine="by_coords")
    if "depth" in ds.dims:
        ds = ds.isel(depth=0)

    return FieldSet.from_data(
        {"U": ds["uo"].values, "V": ds["vo"].values},
        {
            "lon": ds["longitude"].values,
            "lat": ds["latitude"].values,
            "time": ds["time"].values,
        },
        mesh="spherical",
        # A particle that leaves the domain has no velocity to follow. Parcels
        # deletes it, which is the same convention as the NA in derivoce.
        allow_time_extrapolation=False,
    )


def seed(n: int) -> tuple[np.ndarray, np.ndarray]:
    f = np.linspace(0, 1, n)
    lon = RELEASE[0][0] + f * (RELEASE[1][0] - RELEASE[0][0])
    lat = RELEASE[0][1] + f * (RELEASE[1][1] - RELEASE[0][1])
    return lon, lat


def reached(traj: xr.Dataset, region: dict) -> np.ndarray:
    """True for each particle that entered the region at any point."""
    inside = (
        (traj["lon"] >= region["lon"][0]) & (traj["lon"] <= region["lon"][1])
        & (traj["lat"] >= region["lat"][0]) & (traj["lat"] <= region["lat"][1])
    )
    return inside.any(dim="obs").values


def run_release(fieldset: FieldSet, when: datetime, n: int,
                outdir: Path) -> dict:
    lon, lat = seed(n)
    pset = ParticleSet(fieldset=fieldset, pclass=JITParticle,
                       lon=lon, lat=lat, time=np.datetime64(when))

    zarr_path = outdir / f"traj_{when:%Y%m}.zarr"
    out = pset.ParticleFile(name=str(zarr_path), outputdt=timedelta(days=5))
    pset.execute(AdvectionRK4, runtime=timedelta(days=TRACK_DAYS),
                 dt=timedelta(hours=6), output_file=out)

    traj = xr.open_zarr(zarr_path)
    lab = reached(traj, LABRADOR).sum()
    sco = reached(traj, SCOTIAN).sum()
    return dict(YEAR=when.year, MONTH=when.month, released=n,
                labrador=int(lab), scotian=int(sco),
                index=float(lab - sco) / n)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--start", required=True, help="first release, YYYY-MM-DD")
    p.add_argument("--end", required=True, help="last release, YYYY-MM-DD")
    p.add_argument("--particles", type=int, default=200)
    p.add_argument("--every-days", type=int, default=30,
                   help="release interval; the paper used 7")
    args = p.parse_args()

    OUT.mkdir(exist_ok=True)
    start = datetime.fromisoformat(args.start)
    end = datetime.fromisoformat(args.end)

    # Trajectories run three years past the last release, so the fieldset has to
    # cover that too.
    years = list(range(start.year, end.year + 4))
    print(f"loading velocities for {years[0]}-{years[-1]}")
    fieldset = build_fieldset(years)

    rows, when = [], start
    while when <= end:
        print(f"  release {when:%Y-%m} ...", flush=True)
        rows.append(run_release(fieldset, when, args.particles, OUT))
        print(f"    labrador {rows[-1]['labrador']}  "
              f"scotian {rows[-1]['scotian']}  index {rows[-1]['index']:.3f}")
        when += timedelta(days=args.every_days)

    df = pd.DataFrame(rows)
    csv = OUT / "lcr_extended.csv"
    df.to_csv(csv, index=False)
    print(f"\nwrote {csv}  ({len(df)} releases)")


if __name__ == "__main__":
    main()
