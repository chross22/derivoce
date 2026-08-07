"""Fetch daily GLORYS surface velocities for the LCR reproduction.

One file per year, because the full record is ~4 GB and a single request that
size is fragile. Already-fetched years are skipped, so this is resumable.
"""
import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

DATASET = "cmems_mod_glo_phy_my_0.083deg_P1D-m"

# Spans the release line on the Labrador Shelf and the Scotian Shelf, with room
# for particles to wander. Trimming this is the easiest way to lose trajectories.
#
# The eastern edge was at 45W and that was too close. Particles crossing it did
# so at a median of 130 days, while arrivals on the Scotian Shelf take 410, and
# 72% of those losses left at or below 48N -- slope water that could still have
# turned southwest, not North Atlantic Current export. Losing them truncated the
# Scotian term while leaving the Labrador term complete, which biases the
# difference rather than just thinning it.
BOX = dict(lon=(-70.0, -30.0), lat=(36.0, 58.0))

DATA = Path(__file__).parent / "data"

# GLORYS levels are fixed and unevenly spaced, and a window falling between two
# of them returns nothing. Ask for a narrow band around the level itself rather
# than a round number: 0.494 is the surface, 92.326 is the one nearest 100 m.
SURFACE = 0.494
DEPTH_TOLERANCE = 0.02


def velocity_file(year: int, depth: float) -> Path:
    """Depth is part of the name, so runs at different levels do not collide."""
    return DATA / f"glorys_{year}_{depth:g}m.nc"


def fetch_year(year: int, depth: float, overwrite: bool = False) -> Path:
    DATA.mkdir(exist_ok=True)
    out = velocity_file(year, depth)
    if out.exists() and out.stat().st_size > 0 and not overwrite:
        print(f"  {year}: already have it, skipping")
        return out

    cmd = [
        "copernicusmarine", "subset",
        "--dataset-id", DATASET,
        "--variable", "uo", "--variable", "vo",
        "--minimum-longitude", str(BOX["lon"][0]),
        "--maximum-longitude", str(BOX["lon"][1]),
        "--minimum-latitude", str(BOX["lat"][0]),
        "--maximum-latitude", str(BOX["lat"][1]),
        "--minimum-depth", str(depth * (1 - DEPTH_TOLERANCE)),
        "--maximum-depth", str(depth * (1 + DEPTH_TOLERANCE)),
        "--start-datetime", f"{year}-01-01",
        "--end-datetime", f"{year}-12-31",
        "-o", str(DATA), "--output-filename", out.name, "--overwrite",
    ]
    print(f"  {year}: fetching ...", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr[-2000:], file=sys.stderr)
        raise SystemExit(f"fetch failed for {year}")
    print(f"  {year}: {out.stat().st_size / 1e6:.0f} MB")
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--start", required=True, help="YYYY-MM-DD")
    p.add_argument("--end", required=True, help="YYYY-MM-DD")
    p.add_argument("--depth", type=float, default=SURFACE,
                   help="metres; must be a GLORYS level. 0.494 is the surface, "
                        "92.326 is the level nearest 100 m")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--jobs", type=int, default=4,
                   help="years to fetch at once; each is a separate request")
    args = p.parse_args()

    years = list(range(int(args.start[:4]), int(args.end[:4]) + 1))
    print(f"fetching {len(years)} year(s) at {args.depth} m, {args.jobs} at a time")

    # The years are independent requests, so they overlap. Failures are collected
    # rather than raised, so one bad year does not discard the others' downloads.
    failed = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(fetch_year, y, args.depth, args.overwrite): y
                   for y in years}
        for future in as_completed(futures):
            year = futures[future]
            try:
                future.result()
            except SystemExit as err:
                failed.append(year)
                print(f"  {year}: FAILED ({err})", file=sys.stderr)

    if failed:
        raise SystemExit(f"failed years: {sorted(failed)}. Rerun to resume.")
    print("done")


if __name__ == "__main__":
    main()
