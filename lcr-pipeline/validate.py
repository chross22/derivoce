"""Compare the recomputed index against the published series.

Nothing downstream should be trusted until this shows agreement. A
reimplementation with a different integrator, release cadence, and arrival
criterion can easily produce a plausible series that is not the same quantity.

The headline number is the correlation of DESEASONALIZED anomalies, not of the
raw monthly values. Both series carry the annual cycle of the Labrador Current,
and two series that share a seasonal cycle correlate highly whether or not they
measure the same thing -- a raw correlation can clear 0.7 on seasonality alone
and say nothing about retroflection. Significance is reported against an
effective sample size, because consecutive months are not independent draws and
the usual t-test on r would overstate confidence severalfold.

    python validate.py --depth 92.326
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).parent
PUBLISHED = HERE / "data" / "lcr_published.csv"

# Correlations of a short series are unstable; below this many overlapping
# months the verdict is withheld rather than reported with false precision.
MIN_MONTHS = 24
# Deseasonalizing needs enough years for a monthly climatology to mean anything.
MIN_YEARS_FOR_CLIMATOLOGY = 5


def load_pair(published: Path, recomputed: Path) -> pd.DataFrame:
    if not Path(published).exists():
        raise SystemExit(f"missing {published}")
    if not Path(recomputed).exists():
        # Single runs are named for the span they cover; the unspanned name is
        # the merged record that run_record.py writes.
        near = sorted(Path(recomputed).parent.glob("lcr_extended*.csv"))
        hint = ("\n  available:\n    " + "\n    ".join(f.name for f in near)
                if near else "\n  nothing in output/ yet; run track.py first")
        raise SystemExit(f"missing {recomputed}{hint}")

    pub = pd.read_csv(published)
    new = pd.read_csv(recomputed)
    pub.columns = [c.upper() for c in pub.columns]
    new.columns = [c.upper() for c in new.columns]

    merged = pub.merge(new, on=["YEAR", "MONTH"], how="inner")
    if merged.empty:
        raise SystemExit("no overlapping months; nothing to compare")

    dup = merged.duplicated(["YEAR", "MONTH"]).sum()
    if dup:
        raise SystemExit(
            f"{dup} duplicate YEAR/MONTH rows after merging. The recomputed "
            "series must be one row per calendar month -- see to_monthly().")

    merged = merged.sort_values(["YEAR", "MONTH"]).reset_index(drop=True)
    ok = np.isfinite(merged["LCR"]) & np.isfinite(merged["INDEX"])
    return merged[ok].reset_index(drop=True)


def deseasonalize(df: pd.DataFrame, col: str) -> np.ndarray:
    """Subtract the month-of-year mean, leaving anomalies.

    Computed on the overlap itself rather than on a fixed reference period, so
    both series lose the same climatology and neither is advantaged.
    """
    return (df[col] - df.groupby("MONTH")[col].transform("mean")).to_numpy()


def detrend(x: np.ndarray) -> np.ndarray:
    t = np.arange(len(x), dtype=float)
    slope, intercept = np.polyfit(t, x, 1)
    return x - (slope * t + intercept)


def lag1(x: np.ndarray) -> float:
    if len(x) < 3:
        return 0.0
    r = np.corrcoef(x[:-1], x[1:])[0, 1]
    return 0.0 if not np.isfinite(r) else float(r)


def effective_n(a: np.ndarray, b: np.ndarray) -> float:
    """Bretherton et al. (1999): months are not independent draws.

    Two series each autocorrelated at lag one carry far less information than
    their length suggests, so the degrees of freedom behind any p-value have to
    be discounted accordingly.
    """
    r1, r2 = lag1(a), lag1(b)
    factor = (1 - r1 * r2) / (1 + r1 * r2)
    return max(3.0, len(a) * min(1.0, max(factor, 1e-3)))


def p_from_r(r: float, n_eff: float) -> float:
    if not np.isfinite(r) or n_eff <= 2 or abs(r) >= 1:
        return float("nan")
    t = r * np.sqrt((n_eff - 2) / (1 - r ** 2))
    return float(2 * stats.t.sf(abs(t), n_eff - 2))


def block_bootstrap_ci(a: np.ndarray, b: np.ndarray, block: int = 12,
                       draws: int = 2000, seed: int = 0) -> tuple:
    """Confidence interval that respects autocorrelation.

    Resampling single months would break the serial dependence and give an
    interval far too narrow, so contiguous blocks are resampled instead.
    """
    n = len(a)
    if n < 2 * block:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    starts_max = n - block
    out = np.empty(draws)
    for i in range(draws):
        idx = np.concatenate([np.arange(s, s + block)
                              for s in rng.integers(0, starts_max + 1,
                                                    n // block + 1)])[:n]
        with np.errstate(invalid="ignore"):
            out[i] = np.corrcoef(a[idx], b[idx])[0, 1]
    out = out[np.isfinite(out)]
    if out.size == 0:
        return (float("nan"), float("nan"))
    return tuple(np.percentile(out, [2.5, 97.5]))


def report_pair(name: str, a: np.ndarray, b: np.ndarray) -> float:
    r = float(np.corrcoef(a, b)[0, 1])
    rho = float(stats.spearmanr(a, b).statistic)
    n_eff = effective_n(a, b)
    p = p_from_r(r, n_eff)
    lo, hi = block_bootstrap_ci(a, b)
    print(f"  {name:<16} r {r:+.3f}   rho {rho:+.3f}   "
          f"n_eff {n_eff:5.1f}   p {p:.3f}   95% CI [{lo:+.2f}, {hi:+.2f}]")
    return r


def lagged(a: np.ndarray, b: np.ndarray, max_lag: int = 12) -> tuple:
    """Scan for a systematic offset.

    Counting particles into a box and counting them across a section need not
    peak in the same month; a displaced peak is a timing difference, not a
    disagreement about the physics.

    A positive lag means the first series leads the second: a value in `a` at
    month i reappears in `b` at month i + lag.
    """
    best_lag, best_r, rows = 0, -2.0, []
    for lag in range(-max_lag, max_lag + 1):
        if lag < 0:
            x, y = a[-lag:], b[:len(b) + lag]
        elif lag > 0:
            x, y = a[:len(a) - lag], b[lag:]
        else:
            x, y = a, b
        if len(x) < MIN_MONTHS // 2:
            continue
        r = float(np.corrcoef(x, y)[0, 1])
        rows.append((lag, r))
        if r > best_r:
            best_lag, best_r = lag, r
    return best_lag, best_r, rows


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--depth", type=float, default=0.494)
    p.add_argument("--published", default=str(PUBLISHED))
    p.add_argument("--recomputed", default=None)
    args = p.parse_args()

    recomputed = (Path(args.recomputed) if args.recomputed
                  else HERE / "output" / f"lcr_extended_{args.depth:g}m.csv")
    df = load_pair(Path(args.published), recomputed)

    n = len(df)
    years = df["YEAR"].nunique()
    span = f"{df['YEAR'].iloc[0]}-{df['MONTH'].iloc[0]:02d} to " \
           f"{df['YEAR'].iloc[-1]}-{df['MONTH'].iloc[-1]:02d}"
    print(f"{recomputed.name}  vs  published")
    print(f"overlap: {n} months over {years} year(s), {span}\n")

    print(f"published   mean {df['LCR'].mean():+.4f}  sd {df['LCR'].std():.4f}")
    print(f"recomputed  mean {df['INDEX'].mean():+.4f}  sd {df['INDEX'].std():.4f}")
    print("(scales are not expected to match: this counts particles entering "
          "boxes,\n theirs counts crossings of sections)\n")

    raw_a, raw_b = df["LCR"].to_numpy(), df["INDEX"].to_numpy()

    print("CORRELATIONS")
    report_pair("raw", raw_a, raw_b)

    enough = years >= MIN_YEARS_FOR_CLIMATOLOGY
    if enough:
        anom_a = deseasonalize(df, "LCR")
        anom_b = deseasonalize(df, "INDEX")
        r_headline = report_pair("deseasonalized", anom_a, anom_b)
        report_pair("detrended", detrend(raw_a), detrend(raw_b))
        report_pair("both", detrend(anom_a), detrend(anom_b))
    else:
        r_headline = None
        print(f"  deseasonalized   withheld: {years} year(s) of overlap is too "
              f"few for a\n                   monthly climatology (need "
              f"{MIN_YEARS_FOR_CLIMATOLOGY}). The raw number above is\n"
              f"                   NOT evidence the series agree -- a shared "
              f"annual cycle\n                   alone can produce it.")

    if n >= MIN_MONTHS:
        best_lag, best_r, rows = lagged(raw_a, raw_b)
        print(f"\nLAGS  (recomputed shifted against published)")
        for lag, r in rows:
            if lag % 3 == 0 or lag == best_lag:
                mark = "  <- best" if lag == best_lag else ""
                print(f"  {lag:+3d} months  r {r:+.3f}{mark}")
        if abs(best_lag) > 2:
            print(f"  peak at {best_lag:+d} months suggests a systematic timing "
                  "offset, not disagreement")

    print("\nVERDICT")
    if n < MIN_MONTHS:
        print(f"  withheld. {n} overlapping months is too few to correlate "
              f"(need {MIN_MONTHS}).\n  Extend the release period and rerun.")
        return
    if r_headline is None:
        print(f"  withheld. Needs {MIN_YEARS_FOR_CLIMATOLOGY}+ years of overlap "
              "to separate the seasonal\n  cycle from the signal.")
        return

    if r_headline > 0.7:
        print("  tracks the published series; the extension is worth trusting")
    elif r_headline > 0.4:
        print("  partly tracks; revisit the arrival regions before extending")
    else:
        print("  does not track. Do not extend. Run diagnose.py first: check "
              "particle\n  survival and the exit boundary before changing the "
              "boxes, since a\n  truncated ensemble cannot agree with anything")


if __name__ == "__main__":
    main()
