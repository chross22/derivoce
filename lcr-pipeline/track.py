"""Track particles from the Labrador Shelf and count where they arrive.

Seeds a line across the Labrador Shelf at intervals, advects each release for
three years through daily GLORYS velocities, and records how many reach the
Labrador Shelf and how many reach the Scotian Shelf.

The index is the difference between those counts, scaled by the number released.
High retroflection means fewer particles carry on southwest, so the difference
grows.

Every release goes into one ParticleSet with its own start time, and the whole
set is advected in a single compiled pass. Releases are independent, so running
them one at a time would recompile the kernel and re-walk the same velocity
fields for each; staggering the start times inside one set does not.
"""
import argparse
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from parcels import (AdvectionRK4, FieldSet, JITParticle, ParticleSet,
                     StatusCode, Variable)

from fetch import SURFACE, velocity_file

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


class LCRParticle(JITParticle):
    """Carries which release it belongs to, so one set holds them all."""

    release = Variable("release", dtype=np.int32, initial=0, to_write="once")


def DeleteOnError(particle, fieldset, time):  # noqa: N802 (parcels kernel naming)
    """Retire a particle that leaves the domain or outruns the velocities.

    Without this the first particle to drift out of the box aborts the whole
    run. Leaving is a normal outcome here, and matches the NA convention in
    derivoce: the trajectory simply stops.
    """
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle.delete()
    elif particle.state == StatusCode.ErrorThroughSurface:
        particle.delete()
    elif particle.state == StatusCode.ErrorTimeExtrapolation:
        particle.delete()


def build_fieldset(years: list[int], depth: float = SURFACE) -> FieldSet:
    """One fieldset spanning every year needed, so trajectories cross year ends."""
    files = [velocity_file(y, depth) for y in years]
    missing = [f.name for f in files if not f.exists()]
    if missing:
        raise SystemExit(
            f"missing velocity files: {missing}. Run "
            f"fetch.py --depth {depth:g} first.")

    ds = xr.open_mfdataset([str(f) for f in files], combine="by_coords")
    if "depth" in ds.dims:
        if ds.sizes["depth"] > 1:
            raise SystemExit(
                f"expected one depth level, got {ds['depth'].values}. "
                "Narrow DEPTH_TOLERANCE in fetch.py and refetch.")
        ds = ds.isel(depth=0)
    print(f"  level {float(ds['depth'].values):.3f} m, "
          f"{ds.sizes['time']} daily fields")
    report_mask(ds)

    return FieldSet.from_data(
        {"U": ds["uo"].values, "V": ds["vo"].values},
        {
            "lon": ds["longitude"].values,
            "lat": ds["latitude"].values,
            "time": ds["time"].values,
        },
        mesh="spherical",
        # Running past the last velocity field would silently freeze the flow,
        # so parcels raises instead and DeleteOnError retires the particle.
        allow_time_extrapolation=False,
    )


def masked_fraction(ds: xr.Dataset, region: dict) -> float:
    """Share of a region that is land at this depth.

    GLORYS masks each level against the bathymetry, so a level deeper than the
    shelf turns the shelf into land. Particles cannot cross it, and an arrival
    box sitting on it can never fill.
    """
    sub = ds["uo"].isel(time=0).sel(
        longitude=slice(*region["lon"]), latitude=slice(*region["lat"]))
    return float(np.isnan(sub.values).mean())


def report_mask(ds: xr.Dataset, limit: float = 0.25) -> None:
    """Say up front whether the arrival boxes are reachable at this depth.

    A run at a level that buries the shelf produces a clean-looking series of
    near-zero counts, which is indistinguishable from a real absence of
    transport unless the mask is checked.
    """
    for name, region in (("LABRADOR", LABRADOR), ("SCOTIAN", SCOTIAN)):
        frac = masked_fraction(ds, region)
        flag = ("  <- much of this box is land at this depth; a low arrival "
                "count here means little") if frac > limit else ""
        print(f"  {name} box {100 * frac:.1f}% masked{flag}")


def check_coverage(fieldset: FieldSet, dates: list[datetime]) -> None:
    """Fail before advecting if the velocities stop short of the last release.

    A trajectory that outruns the fields is retired by DeleteOnError, which is
    right for a particle leaving the domain but would quietly truncate every
    late release if the fields simply end early.
    """
    last_needed = dates[-1] + timedelta(days=TRACK_DAYS)
    have = fieldset.U.grid.time_full
    end = fieldset.U.grid.time_origin.fulltime(have[-1])
    end = pd.Timestamp(end).to_pydatetime()
    if end < last_needed:
        short = (last_needed - end).days
        raise SystemExit(
            f"velocities end {end:%Y-%m-%d}, but the last release needs them to "
            f"{last_needed:%Y-%m-%d} ({short} days short). Fetch through "
            f"{last_needed.year}."
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


def release_dates(start: datetime, end: datetime, every_days: int) -> list[datetime]:
    dates, when = [], start
    while when <= end:
        dates.append(when)
        when += timedelta(days=every_days)
    return dates


def run_all(fieldset: FieldSet, dates: list[datetime], n: int,
            outdir: Path, tag: str = "") -> Path:
    """Advect every release together, and write one trajectory file."""
    lon, lat = seed(n)
    all_lon = np.tile(lon, len(dates))
    all_lat = np.tile(lat, len(dates))
    all_time = np.repeat([np.datetime64(d) for d in dates], n)
    all_release = np.repeat(np.arange(len(dates), dtype=np.int32), n)

    pset = ParticleSet(fieldset=fieldset, pclass=LCRParticle,
                       lon=all_lon, lat=all_lat, time=all_time,
                       release=all_release)

    # The last release still needs its full three years, so the run extends
    # past the final start time rather than three years from the first.
    runtime = (dates[-1] - dates[0]) + timedelta(days=TRACK_DAYS)

    zarr_path = outdir / f"traj_{dates[0]:%Y%m}_{dates[-1]:%Y%m}{tag}.zarr"
    out = pset.ParticleFile(name=str(zarr_path), outputdt=timedelta(days=5))
    print(f"  advecting {len(all_lon)} particles over {runtime.days} days", flush=True)
    pset.execute([AdvectionRK4, DeleteOnError], runtime=runtime,
                 dt=timedelta(hours=6), output_file=out)
    return zarr_path


def tally_releases(traj: xr.Dataset, dates: list[datetime],
                   labrador: dict = LABRADOR,
                   scotian: dict = SCOTIAN) -> pd.DataFrame:
    """Count arrivals for each release.

    Takes an open dataset and the two regions as arguments so the arrival boxes
    can be swept without advecting again: the trajectories do not depend on
    where the boxes are, only the counting does.
    """
    in_lab = reached(traj, labrador)
    in_sco = reached(traj, scotian)
    which = np.asarray(traj["release"].values).reshape(-1)

    rows = []
    for i, when in enumerate(dates):
        mine = which == i
        rows.append(dict(DATE=when, YEAR=when.year, MONTH=when.month,
                         released=int(mine.sum()),
                         labrador=int(in_lab[mine].sum()),
                         scotian=int(in_sco[mine].sum())))
    return pd.DataFrame(rows)


def to_monthly(per_release: pd.DataFrame) -> pd.DataFrame:
    """Collapse releases onto calendar months.

    Releases are on a fixed day interval, so they drift against the calendar:
    at weekly cadence most months hold four and some five, and at 30 days a
    month can hold two or none. The published series is monthly, and
    validate.py merges on YEAR/MONTH, so releases have to be pooled per month
    rather than assumed one-to-one -- otherwise a month with two releases
    yields a duplicate key and a month with none silently drops out.

    Counts are pooled before dividing, so a month with more releases is not
    given more weight than its particles earn.
    """
    grouped = (per_release.groupby(["YEAR", "MONTH"], as_index=False)
               .agg(releases=("released", "size"), released=("released", "sum"),
                    labrador=("labrador", "sum"), scotian=("scotian", "sum")))
    grouped["index"] = ((grouped["labrador"] - grouped["scotian"])
                        / grouped["released"].where(grouped["released"] > 0))
    return grouped.sort_values(["YEAR", "MONTH"]).reset_index(drop=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--start", required=True, help="first release, YYYY-MM-DD")
    p.add_argument("--end", required=True, help="last release, YYYY-MM-DD")
    p.add_argument("--particles", type=int, default=200)
    p.add_argument("--every-days", type=int, default=7,
                   help="release interval in days; the paper used 7")
    p.add_argument("--depth", type=float, default=SURFACE,
                   help="GLORYS level to advect on; must already be fetched")
    args = p.parse_args()

    OUT.mkdir(exist_ok=True)
    start = datetime.fromisoformat(args.start)
    end = datetime.fromisoformat(args.end)

    # Trajectories run three years past the last release, so the fieldset has to
    # cover that too.
    years = list(range(start.year, end.year + 4))
    print(f"loading velocities for {years[0]}-{years[-1]} at {args.depth:g} m")
    fieldset = build_fieldset(years, args.depth)

    dates = release_dates(start, end, args.every_days)
    check_coverage(fieldset, dates)
    print(f"{len(dates)} releases of {args.particles} particles, "
          f"every {args.every_days} days")
    tag = f"_{args.depth:g}m"
    zarr_path = run_all(fieldset, dates, args.particles, OUT, tag)

    per_release = tally_releases(xr.open_zarr(zarr_path), dates)
    df = to_monthly(per_release)
    for row in df.itertuples():
        print(f"  {row.YEAR}-{row.MONTH:02d}  {row.releases} release(s)  "
              f"labrador {row.labrador}  scotian {row.scotian}  "
              f"index {row.index:.3f}")

    per_release.to_csv(OUT / f"lcr_releases{tag}.csv", index=False)
    csv = OUT / f"lcr_extended{tag}.csv"
    df.to_csv(csv, index=False)
    print(f"\nwrote {csv}  ({len(df)} months from {len(per_release)} releases)")


if __name__ == "__main__":
    main()
