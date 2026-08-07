"""Fetch daily GLORYS surface velocities for the LCR reproduction.

One file per year, because the full record is ~4 GB and a single request that
size is fragile. Already-fetched years are skipped, so this is resumable.
"""
import argparse
import subprocess
import sys
from pathlib import Path

DATASET = "cmems_mod_glo_phy_my_0.083deg_P1D-m"

# Spans the release line on the Labrador Shelf and the Scotian Shelf, with room
# for particles to wander. Trimming this is the easiest way to lose trajectories.
BOX = dict(lon=(-70.0, -45.0), lat=(40.0, 58.0))

DATA = Path(__file__).parent / "data"


def fetch_year(year: int, depth: float, overwrite: bool = False) -> Path:
    DATA.mkdir(exist_ok=True)
    out = DATA / f"glorys_{year}.nc"
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
        "--minimum-depth", str(depth), "--maximum-depth", str(depth + 1),
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
    p.add_argument("--depth", type=float, default=0.49,
                   help="metres; 0.49 is the surface level")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args()

    years = range(int(args.start[:4]), int(args.end[:4]) + 1)
    print(f"fetching {len(list(years))} year(s) at {args.depth} m")
    for year in years:
        fetch_year(year, args.depth, args.overwrite)
    print("done")


if __name__ == "__main__":
    main()
