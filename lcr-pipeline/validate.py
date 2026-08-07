"""Compare the recomputed index against the published series.

Nothing downstream should be trusted until this shows agreement. A
reimplementation with a different integrator, release cadence, and arrival
criterion can easily produce a plausible series that is not the same quantity.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
PUBLISHED = HERE / "data" / "lcr_published.csv"
RECOMPUTED = HERE / "output" / "lcr_extended.csv"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--published", default=str(PUBLISHED),
                   help="CSV with YEAR, MONTH, LCR. Export from R with: "
                        "write.csv(datamatch::fetch_climate_index('LCR'), ...)")
    p.add_argument("--recomputed", default=str(RECOMPUTED))
    args = p.parse_args()

    for path in (args.published, args.recomputed):
        if not Path(path).exists():
            raise SystemExit(f"missing {path}")

    pub = pd.read_csv(args.published)
    new = pd.read_csv(args.recomputed)
    pub.columns = [c.upper() for c in pub.columns]
    new.columns = [c.upper() for c in new.columns]

    merged = pub.merge(new, on=["YEAR", "MONTH"], suffixes=("_pub", "_new"))
    if merged.empty:
        raise SystemExit("no overlapping months; nothing to compare")

    a = merged["LCR"].to_numpy()
    b = merged["INDEX"].to_numpy()
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]

    print(f"overlapping months: {len(a)}")
    if len(a) < 3:
        raise SystemExit("too few overlapping months to say anything")

    r = np.corrcoef(a, b)[0, 1]
    print(f"correlation: {r:.3f}")
    print(f"published  mean {a.mean():.4f}  sd {a.std():.4f}")
    print(f"recomputed mean {b.mean():.4f}  sd {b.std():.4f}")

    # The scales will differ: this counts particles reaching boxes, theirs counts
    # crossings of hydrographic sections. Correlation is the question, not the
    # absolute value, so the sensible check is whether the two move together.
    if r > 0.7:
        verdict = "tracks the published series; the extension is worth trusting"
    elif r > 0.4:
        verdict = "partly tracks; revisit the arrival regions before extending"
    else:
        verdict = "does not track. Do not extend. Check the arrival regions, " \
                  "release cadence, and depth before anything else"
    print(f"\nverdict: {verdict}")


if __name__ == "__main__":
    main()
