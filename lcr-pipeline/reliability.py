"""How much of the index is signal, and how much is counting noise?

A failure to correlate with the published series has two very different causes.
Either the method measures the wrong thing, or it measures the right thing too
noisily to tell. This separates them without any reference to the published
series at all: split each month's particles into halves, compute the index
independently on each, and see whether the index agrees with itself.

The paper releases 966 particles a week; these runs release 100, so counting
noise is a real candidate. Split-half reliability caps how well the index can
correlate with anything -- the ceiling is sqrt(reliability). If the ceiling is
high and the observed correlation is low, more particles will not help.

    python reliability.py --depth 47.374
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from sweep_boxes import inside, load_blocks, monthly
from track import LABRADOR, SCOTIAN
from validate import deseasonalize

HERE = Path(__file__).parent


def split_half(dates, lab, sco, seed: int) -> float:
    """Correlation between two independent halves of the same particles."""
    rng = np.random.default_rng(seed)
    pick = rng.random(len(dates)) < 0.5

    a = monthly(dates[pick], lab[pick], sco[pick])
    b = monthly(dates[~pick], lab[~pick], sco[~pick])
    m = a.merge(b, on=["YEAR", "MONTH"], suffixes=("_a", "_b"))
    if len(m) < 24:
        return np.nan

    x = deseasonalize(m.rename(columns={"INDEX_a": "LCR"}), "LCR")
    y = deseasonalize(m.rename(columns={"INDEX_b": "LCR"}), "LCR")
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--depth", type=float, default=47.374)
    p.add_argument("--draws", type=int, default=200)
    args = p.parse_args()

    print("loading trajectories")
    lon, lat, dates = load_blocks(args.depth)
    lab = inside(lon, lat, LABRADOR)
    sco = inside(lon, lat, SCOTIAN)
    print(f"  {len(dates)} particles, {dates.min():%Y-%m} to {dates.max():%Y-%m}\n")

    rs = np.array([split_half(dates, lab, sco, s) for s in range(args.draws)])
    rs = rs[np.isfinite(rs)]
    if rs.size == 0:
        raise SystemExit("could not form split halves")

    half = float(np.median(rs))
    # Spearman-Brown: each half has half the particles, so the full-sample
    # reliability is higher than the agreement between halves.
    full = 2 * half / (1 + half) if half > -1 else np.nan
    ceiling = np.sqrt(full) if full > 0 else 0.0

    print(f"split-half agreement over {rs.size} random splits")
    print(f"  median r        {half:+.3f}   (IQR {np.percentile(rs, 25):+.3f} "
          f"to {np.percentile(rs, 75):+.3f})")
    print(f"  full-sample reliability (Spearman-Brown)  {full:+.3f}")
    print(f"  ceiling on any correlation                 {ceiling:.3f}")

    print("\ninterpretation")
    if full < 0.3:
        print("  The index barely agrees with itself. It is dominated by counting")
        print("  noise at this particle count, and no correlation with the")
        print("  published series is possible until that is fixed. Raise")
        print("  --particles and rerun before concluding anything about method.")
    elif ceiling < 0.7:
        print(f"  The index is only moderately reproducible. Even a perfect")
        print(f"  method could not exceed r = {ceiling:.2f} here, short of the 0.7")
        print("  target, so more particles are needed regardless.")
    else:
        print(f"  The index agrees with itself well: a correct method could")
        print(f"  reach r = {ceiling:.2f}. Counting noise is not the limitation,")
        print("  so a low correlation against the published series is a")
        print("  statement about the method, not the sample size.")


if __name__ == "__main__":
    main()
