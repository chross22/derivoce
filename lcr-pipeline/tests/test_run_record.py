"""Tests for the block splitting used by the parallel record run.

An off-by-one here does not crash: it runs a block against velocities that stop
early, and DeleteOnError quietly retires every late trajectory. The result is a
plausible series computed from truncated tracks, so the arithmetic is checked
rather than trusted.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import run_record
import track


class TestBlockYears:
    def test_a_block_covers_the_full_tracking_window(self):
        """Releases through December of Y are tracked TRACK_DAYS past it."""
        years = run_record.block_years(1995)
        last_release = 1995 + 1  # a December release runs into the next years
        needed_through = 1995 + track.TRACK_DAYS // 365
        assert years[0] == 1995
        assert years[-1] >= needed_through
        assert last_release in years

    def test_matches_what_track_py_asks_for(self):
        """track.py builds its fieldset from range(start.year, end.year + 4).

        The two must agree, or run_record will report files as present that
        track.py then rejects, or vice versa.
        """
        for year in (1993, 1995, 2014):
            assert run_record.block_years(year) == list(range(year, year + 4))

    def test_blocks_are_contiguous_and_overlapping(self):
        a = run_record.block_years(1995)
        b = run_record.block_years(1996)
        assert b[0] == a[0] + 1
        assert set(a) & set(b), "consecutive blocks share velocity years"


class TestMissingFiles:
    def test_reports_every_year_a_block_needs(self, tmp_path, monkeypatch):
        monkeypatch.setattr(run_record, "velocity_file",
                            lambda y, d: tmp_path / f"glorys_{y}_{d:g}m.nc")
        gaps = run_record.missing_files(range(1995, 1996), 47.374)
        assert len(gaps) == 4
        assert "glorys_1998_47.374m.nc" in gaps

    def test_silent_when_everything_is_present(self, tmp_path, monkeypatch):
        monkeypatch.setattr(run_record, "velocity_file",
                            lambda y, d: tmp_path / f"glorys_{y}_{d:g}m.nc")
        for y in range(1995, 1999):
            (tmp_path / f"glorys_{y}_47.374m.nc").touch()
        assert run_record.missing_files(range(1995, 1996), 47.374) == []

    def test_shared_years_are_not_double_reported(self, tmp_path, monkeypatch):
        """Consecutive blocks overlap; a gap is one file, not one per block."""
        monkeypatch.setattr(run_record, "velocity_file",
                            lambda y, d: tmp_path / f"glorys_{y}_{d:g}m.nc")
        gaps = run_record.missing_files(range(1995, 1998), 47.374)
        assert len(gaps) == len(set(gaps))
        assert len(gaps) == 6  # 1995..2000
