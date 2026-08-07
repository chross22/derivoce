"""Run the whole record in parallel blocks, and merge the result.

Releases in year Y need velocities for Y through Y+3 and nothing else, so each
year is an independent job. Running them as separate processes keeps every
fieldset small -- one four-year window rather than the whole record -- and uses
the cores that a single sequential run leaves idle.

Separate processes rather than threads or a pool: parcels compiles a kernel per
ParticleSet and writes C to a scratch directory, and independent interpreters
cannot collide over it.

    python run_record.py --start 1993 --end 2014 --depth 47.374 --jobs 4

Memory is the limit on --jobs, not cores. Each worker holds four years of daily
velocities: roughly 1.5 GB at the current box, so four workers is about 6 GB.
"""
import argparse
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from fetch import SURFACE, velocity_file
from track import TRACK_DAYS, to_monthly

HERE = Path(__file__).parent
OUT = HERE / "output"

# A release in year Y is tracked TRACK_DAYS past the last release in December,
# so the fields have to reach into Y+3.
YEARS_AHEAD = TRACK_DAYS // 365 + 1


def block_years(year: int) -> list[int]:
    return list(range(year, year + YEARS_AHEAD))


def missing_files(years: range, depth: float) -> list[str]:
    needed = sorted({y for year in years for y in block_years(year)})
    return [velocity_file(y, depth).name
            for y in needed if not velocity_file(y, depth).exists()]


def run_block(year: int, depth: float, particles: int, every: int) -> tuple:
    """One year of releases, as its own process."""
    cmd = [sys.executable, str(HERE / "track.py"),
           "--start", f"{year}-01-01", "--end", f"{year}-12-31",
           "--particles", str(particles), "--every-days", str(every),
           "--depth", str(depth)]
    done = subprocess.run(cmd, cwd=str(HERE), capture_output=True, text=True)
    return year, done.returncode, done.stderr[-1500:] if done.returncode else ""


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--start", type=int, required=True, help="first release year")
    p.add_argument("--end", type=int, required=True, help="last release year")
    p.add_argument("--depth", type=float, default=SURFACE)
    p.add_argument("--particles", type=int, default=100)
    p.add_argument("--every-days", type=int, default=7)
    p.add_argument("--jobs", type=int, default=4,
                   help="concurrent blocks; limited by memory, not cores")
    args = p.parse_args()

    years = range(args.start, args.end + 1)
    gaps = missing_files(years, args.depth)
    if gaps:
        raise SystemExit(
            f"missing {len(gaps)} velocity files, first few: {gaps[:4]}.\n"
            f"Releases through {args.end} need fields to {args.end + YEARS_AHEAD - 1}. "
            f"Run fetch.py --depth {args.depth:g}")

    print(f"{len(list(years))} release years, {args.jobs} at a time")
    failed = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_block, y, args.depth, args.particles,
                               args.every_days): y for y in years}
        for future in as_completed(futures):
            year, code, err = future.result()
            if code:
                failed.append(year)
                print(f"  {year}: FAILED\n{err}", file=sys.stderr)
            else:
                print(f"  {year}: done", flush=True)

    if failed:
        raise SystemExit(f"failed years: {sorted(failed)}. Rerun to redo them.")

    tag = f"_{args.depth:g}m"
    parts = sorted(OUT.glob(f"lcr_releases_*{tag}.csv"))
    if not parts:
        raise SystemExit("no per-release output to merge")

    per_release = pd.concat([pd.read_csv(f) for f in parts], ignore_index=True)
    per_release = per_release.drop_duplicates(["DATE"]).sort_values("DATE")
    monthly = to_monthly(per_release)

    csv = OUT / f"lcr_extended{tag}.csv"
    monthly.to_csv(csv, index=False)
    print(f"\nwrote {csv}  ({len(monthly)} months from {len(per_release)} releases)")
    print(f"next: python validate.py --depth {args.depth:g}")


if __name__ == "__main__":
    main()
