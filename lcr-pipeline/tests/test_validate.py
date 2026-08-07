"""Correctness tests for the comparison statistics.

The verdict rests on these, so they are checked against series whose answer is
known by construction: a shared seasonal cycle must not read as agreement, a
known lag must be recovered, and autocorrelation must cost degrees of freedom.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import validate


def monthly_frame(n_years, lcr, index):
    return pd.DataFrame({
        "YEAR": np.repeat(np.arange(1993, 1993 + n_years), 12),
        "MONTH": np.tile(np.arange(1, 13), n_years),
        "LCR": lcr,
        "INDEX": index,
    })


class TestDeseasonalize:
    def test_a_pure_seasonal_cycle_becomes_zero(self):
        n_years = 10
        cycle = np.sin(2 * np.pi * np.arange(1, 13) / 12)
        df = monthly_frame(n_years, np.tile(cycle, n_years), np.tile(cycle, n_years))
        assert validate.deseasonalize(df, "LCR") == pytest.approx(np.zeros(120),
                                                                 abs=1e-12)

    def test_shared_seasonality_does_not_read_as_agreement(self):
        """The trap the headline number exists to avoid.

        Two series built from the same annual cycle plus independent noise
        correlate strongly raw, and not at all once the cycle is removed. If
        this ever fails, the verdict is reporting seasonality as skill.
        """
        rng = np.random.default_rng(0)
        n_years = 20
        cycle = 3 * np.sin(2 * np.pi * np.arange(1, 13) / 12)
        seasonal = np.tile(cycle, n_years)
        a = seasonal + rng.normal(size=240)
        b = seasonal + rng.normal(size=240)
        df = monthly_frame(n_years, a, b)

        raw = np.corrcoef(a, b)[0, 1]
        anom = np.corrcoef(validate.deseasonalize(df, "LCR"),
                           validate.deseasonalize(df, "INDEX"))[0, 1]

        assert raw > 0.7, "expected the shared cycle to inflate the raw number"
        assert abs(anom) < 0.2, "deseasonalizing should remove the false skill"

    def test_genuine_agreement_survives_deseasonalizing(self):
        """The complement: real common signal must not be removed."""
        rng = np.random.default_rng(1)
        n_years = 20
        signal = rng.normal(size=240)
        cycle = np.tile(3 * np.sin(2 * np.pi * np.arange(1, 13) / 12), n_years)
        a = cycle + signal
        b = cycle + signal + 0.3 * rng.normal(size=240)
        df = monthly_frame(n_years, a, b)

        anom = np.corrcoef(validate.deseasonalize(df, "LCR"),
                           validate.deseasonalize(df, "INDEX"))[0, 1]
        assert anom > 0.8


class TestDetrend:
    def test_a_pure_trend_becomes_zero(self):
        x = 5.0 + 2.0 * np.arange(50)
        assert validate.detrend(x) == pytest.approx(np.zeros(50), abs=1e-9)

    def test_result_has_no_trend_left(self):
        """A finite noise sample has some fitted slope; detrending removes it."""
        rng = np.random.default_rng(2)
        x = rng.normal(size=100)
        out = validate.detrend(x)
        slope = np.polyfit(np.arange(100.0), out, 1)[0]
        assert slope == pytest.approx(0.0, abs=1e-12)
        assert out.mean() == pytest.approx(0.0, abs=1e-12)

    def test_structure_is_preserved(self):
        rng = np.random.default_rng(2)
        x = rng.normal(size=100)
        assert np.corrcoef(validate.detrend(x), x)[0, 1] > 0.98


class TestEffectiveN:
    def test_white_noise_keeps_most_of_its_degrees_of_freedom(self):
        rng = np.random.default_rng(3)
        a, b = rng.normal(size=400), rng.normal(size=400)
        assert validate.effective_n(a, b) > 300

    def test_autocorrelation_costs_degrees_of_freedom(self):
        """Persistent series carry less information than their length implies."""
        rng = np.random.default_rng(4)

        def red(n, phi=0.9):
            x = np.zeros(n)
            for i in range(1, n):
                x[i] = phi * x[i - 1] + rng.normal()
            return x

        a, b = red(400), red(400)
        assert validate.effective_n(a, b) < 100

    def test_never_returns_fewer_than_three(self):
        assert validate.effective_n(np.arange(4.0), np.arange(4.0)) >= 3

    def test_discounted_n_widens_the_p_value(self):
        r = 0.35
        assert validate.p_from_r(r, 240) < validate.p_from_r(r, 30)


class TestLagged:
    def test_a_known_offset_is_recovered(self):
        """Positive lag means the first series leads the second.

        a is base advanced by four steps, so a[i] reappears in b at i+4.
        """
        rng = np.random.default_rng(5)
        base = rng.normal(size=200)
        shift = 4
        a = base[shift:]
        b = base[:-shift]
        assert a[0] == pytest.approx(b[shift])

        best_lag, best_r, _ = validate.lagged(a, b, max_lag=12)
        assert best_lag == shift
        assert best_r > 0.95

    def test_the_opposite_offset_gets_the_opposite_sign(self):
        rng = np.random.default_rng(9)
        base = rng.normal(size=200)
        shift = 3
        best_lag, _, _ = validate.lagged(base[:-shift], base[shift:], max_lag=12)
        assert best_lag == -shift

    def test_zero_lag_wins_when_series_are_aligned(self):
        rng = np.random.default_rng(6)
        base = rng.normal(size=200)
        best_lag, _, _ = validate.lagged(base, base + 0.1 * rng.normal(size=200),
                                         max_lag=12)
        assert best_lag == 0


class TestLoadPair:
    def test_duplicate_months_are_rejected(self, tmp_path):
        """A duplicated key would silently mis-pair rows in the merge."""
        pub = tmp_path / "pub.csv"
        new = tmp_path / "new.csv"
        pd.DataFrame({"YEAR": [1995, 1995], "MONTH": [1, 2],
                      "LCR": [0.1, 0.2]}).to_csv(pub, index=False)
        pd.DataFrame({"YEAR": [1995, 1995, 1995], "MONTH": [1, 1, 2],
                      "INDEX": [0.5, 0.6, 0.7]}).to_csv(new, index=False)

        with pytest.raises(SystemExit, match="duplicate"):
            validate.load_pair(pub, new)

    def test_non_overlapping_series_are_rejected(self, tmp_path):
        pub = tmp_path / "pub.csv"
        new = tmp_path / "new.csv"
        pd.DataFrame({"YEAR": [1995], "MONTH": [1], "LCR": [0.1]}).to_csv(
            pub, index=False)
        pd.DataFrame({"YEAR": [2020], "MONTH": [1], "INDEX": [0.5]}).to_csv(
            new, index=False)

        with pytest.raises(SystemExit, match="no overlapping months"):
            validate.load_pair(pub, new)

    def test_rows_come_back_in_date_order(self, tmp_path):
        pub = tmp_path / "pub.csv"
        new = tmp_path / "new.csv"
        frame = pd.DataFrame({"YEAR": [1996, 1995, 1995], "MONTH": [1, 12, 6]})
        frame.assign(LCR=[0.1, 0.2, 0.3]).to_csv(pub, index=False)
        frame.assign(INDEX=[0.4, 0.5, 0.6]).to_csv(new, index=False)

        out = validate.load_pair(pub, new)
        assert out[["YEAR", "MONTH"]].values.tolist() == [[1995, 6], [1995, 12],
                                                          [1996, 1]]


class TestBootstrap:
    def test_interval_brackets_a_strong_correlation(self):
        rng = np.random.default_rng(7)
        a = rng.normal(size=240)
        b = a + 0.3 * rng.normal(size=240)
        lo, hi = validate.block_bootstrap_ci(a, b)
        assert lo < np.corrcoef(a, b)[0, 1] < hi
        assert lo > 0.5

    def test_interval_spans_zero_for_unrelated_series(self):
        rng = np.random.default_rng(8)
        a, b = rng.normal(size=240), rng.normal(size=240)
        lo, hi = validate.block_bootstrap_ci(a, b)
        assert lo < 0 < hi

    def test_short_series_returns_nan_rather_than_a_false_interval(self):
        lo, hi = validate.block_bootstrap_ci(np.arange(10.0), np.arange(10.0))
        assert np.isnan(lo) and np.isnan(hi)
