"""Correctness tests for the tracking machinery.

These do not advect anything. They check the parts that turn trajectories into
a monthly series -- seeding, arrival counting, and calendar aggregation --
against hand-built inputs whose answers are known by inspection. The physics is
checked by validate.py against the published series; this checks the plumbing
that would otherwise make a wrong answer look plausible.

Run from lcr-pipeline/:  pytest tests/
"""
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import track


def trajectories(lon_rows, lat_rows, release):
    """Build a dataset shaped like the zarr parcels writes.

    Rows are particles, columns are observations. NaN stands for a particle
    that was deleted partway, which is what leaving the domain produces.
    """
    return xr.Dataset(
        {
            "lon": (("trajectory", "obs"), np.asarray(lon_rows, dtype=float)),
            "lat": (("trajectory", "obs"), np.asarray(lat_rows, dtype=float)),
            "release": ("trajectory", np.asarray(release, dtype=np.int32)),
        }
    )


class TestSeed:
    def test_endpoints_land_on_the_release_line(self):
        lon, lat = track.seed(10)
        (lon0, lat0), (lon1, lat1) = track.RELEASE
        assert lon[0] == pytest.approx(lon0)
        assert lat[0] == pytest.approx(lat0)
        assert lon[-1] == pytest.approx(lon1)
        assert lat[-1] == pytest.approx(lat1)

    def test_count_is_what_was_asked_for(self):
        lon, lat = track.seed(37)
        assert len(lon) == len(lat) == 37

    def test_spacing_is_even(self):
        lon, _ = track.seed(5)
        steps = np.diff(lon)
        assert steps == pytest.approx(np.full(4, steps[0]))


class TestReleaseDates:
    def test_cadence_and_inclusive_end(self):
        dates = track.release_dates(datetime(1995, 1, 1), datetime(1995, 1, 29), 7)
        assert dates == [datetime(1995, 1, d) for d in (1, 8, 15, 22, 29)]

    def test_end_between_releases_is_not_overshot(self):
        dates = track.release_dates(datetime(1995, 1, 1), datetime(1995, 1, 30), 7)
        assert dates[-1] == datetime(1995, 1, 29)

    def test_single_release_when_start_equals_end(self):
        dates = track.release_dates(datetime(1995, 1, 1), datetime(1995, 1, 1), 7)
        assert dates == [datetime(1995, 1, 1)]

    def test_a_year_of_weekly_releases_covers_every_month(self):
        dates = track.release_dates(datetime(1995, 1, 1), datetime(1995, 12, 31), 7)
        assert sorted({d.month for d in dates}) == list(range(1, 13))


class TestReached:
    region = dict(lon=(-60.0, -50.0), lat=(45.0, 50.0))

    def test_particle_passing_through_counts(self):
        # Second observation is inside; the others are not.
        traj = trajectories([[-70.0, -55.0, -40.0]], [[40.0, 47.0, 40.0]], [0])
        assert track.reached(traj, self.region).tolist() == [True]

    def test_particle_that_never_enters_does_not_count(self):
        traj = trajectories([[-70.0, -65.0]], [[40.0, 41.0]], [0])
        assert track.reached(traj, self.region).tolist() == [False]

    def test_right_longitude_wrong_latitude_does_not_count(self):
        """Both coordinates have to be satisfied at the same observation."""
        traj = trajectories([[-55.0, -70.0]], [[30.0, 47.0]], [0])
        assert track.reached(traj, self.region).tolist() == [False]

    def test_boundary_is_inclusive(self):
        traj = trajectories([[-60.0]], [[45.0]], [0])
        assert track.reached(traj, self.region).tolist() == [True]

    def test_deleted_particle_does_not_count(self):
        """A particle retired mid-run has NaN afterwards, which is not inside."""
        traj = trajectories([[-70.0, np.nan]], [[40.0, np.nan]], [0])
        assert track.reached(traj, self.region).tolist() == [False]

    def test_arrival_before_deletion_still_counts(self):
        traj = trajectories([[-55.0, np.nan]], [[47.0, np.nan]], [0])
        assert track.reached(traj, self.region).tolist() == [True]


class TestTallyReleases:
    labrador = dict(lon=(-60.0, -50.0), lat=(45.0, 50.0))
    scotian = dict(lon=(-70.0, -62.0), lat=(40.0, 44.0))

    def test_counts_are_attributed_to_the_right_release(self):
        # Four particles: two in release 0, two in release 1. One particle of
        # release 0 reaches Labrador; both of release 1 reach Scotian.
        traj = trajectories(
            lon_rows=[[-55.0], [0.0], [-65.0], [-65.0]],
            lat_rows=[[47.0], [0.0], [42.0], [42.0]],
            release=[0, 0, 1, 1],
        )
        dates = [datetime(1995, 1, 1), datetime(1995, 2, 1)]
        out = track.tally_releases(traj, dates, self.labrador, self.scotian)

        assert out["released"].tolist() == [2, 2]
        assert out["labrador"].tolist() == [1, 0]
        assert out["scotian"].tolist() == [0, 2]

    def test_regions_are_arguments_so_boxes_can_be_swept(self):
        """Moving the box changes the count without touching the trajectories."""
        traj = trajectories([[-55.0]], [[47.0]], [0])
        dates = [datetime(1995, 1, 1)]

        hit = track.tally_releases(traj, dates, self.labrador, self.scotian)
        elsewhere = dict(lon=(-30.0, -20.0), lat=(10.0, 20.0))
        miss = track.tally_releases(traj, dates, elsewhere, self.scotian)

        assert hit["labrador"].tolist() == [1]
        assert miss["labrador"].tolist() == [0]


class TestMaskedFraction:
    """A depth level deeper than the shelf turns the shelf into land.

    That produces near-zero arrival counts that look exactly like a real
    absence of transport, so the mask has to be reported, not inferred.
    """

    @staticmethod
    def field(nan_south_of=None):
        lon = np.arange(-70.0, -45.0, 1.0)
        lat = np.arange(40.0, 58.0, 1.0)
        u = np.ones((1, len(lat), len(lon)))
        if nan_south_of is not None:
            u[:, lat < nan_south_of, :] = np.nan
        return xr.Dataset(
            {"uo": (("time", "latitude", "longitude"), u)},
            coords={"time": [0], "latitude": lat, "longitude": lon},
        )

    def test_all_water_is_zero(self):
        region = dict(lon=(-60.0, -55.0), lat=(45.0, 50.0))
        assert track.masked_fraction(self.field(), region) == 0.0

    def test_all_land_is_one(self):
        region = dict(lon=(-60.0, -55.0), lat=(41.0, 44.0))
        assert track.masked_fraction(self.field(nan_south_of=50), region) == 1.0

    def test_partial_mask_is_between(self):
        region = dict(lon=(-60.0, -55.0), lat=(46.0, 53.0))
        frac = track.masked_fraction(self.field(nan_south_of=50), region)
        assert 0.0 < frac < 1.0

    def test_only_the_named_region_is_measured(self):
        """A mask outside the box must not count against it."""
        region = dict(lon=(-60.0, -55.0), lat=(52.0, 56.0))
        assert track.masked_fraction(self.field(nan_south_of=50), region) == 0.0


class TestToMonthly:
    def test_two_releases_in_one_month_pool_into_one_row(self):
        """The bug this guards: 30-day cadence puts two releases in a month."""
        per_release = pd.DataFrame([
            dict(DATE=datetime(1995, 1, 1), YEAR=1995, MONTH=1,
                 released=100, labrador=40, scotian=10),
            dict(DATE=datetime(1995, 1, 31), YEAR=1995, MONTH=1,
                 released=100, labrador=20, scotian=30),
        ])
        out = track.to_monthly(per_release)

        assert len(out) == 1
        assert out.loc[0, "releases"] == 2
        assert out.loc[0, "released"] == 200
        assert out.loc[0, "index"] == pytest.approx((60 - 40) / 200)

    def test_no_duplicate_year_month_keys_survive(self):
        """validate.py merges on YEAR/MONTH, so duplicates would mis-pair."""
        per_release = pd.DataFrame([
            dict(DATE=datetime(1995, 1, d), YEAR=1995, MONTH=1,
                 released=10, labrador=1, scotian=1)
            for d in (1, 8, 15, 22, 29)
        ])
        out = track.to_monthly(per_release)
        assert len(out) == len(out.drop_duplicates(["YEAR", "MONTH"]))

    def test_months_are_pooled_not_averaged(self):
        """A five-release month must not outweigh its particles.

        Averaging the per-release indices would give 0.5; pooling the counts
        gives the fraction of particles that actually arrived.
        """
        per_release = pd.DataFrame([
            dict(DATE=datetime(1995, 1, 1), YEAR=1995, MONTH=1,
                 released=10, labrador=10, scotian=0),
            dict(DATE=datetime(1995, 1, 8), YEAR=1995, MONTH=1,
                 released=90, labrador=0, scotian=0),
        ])
        out = track.to_monthly(per_release)
        assert out.loc[0, "index"] == pytest.approx(10 / 100)

    def test_month_with_no_releases_is_absent_not_zero(self):
        """A gap must not be reported as an index of zero, which is a real value."""
        per_release = pd.DataFrame([
            dict(DATE=datetime(1995, 1, 1), YEAR=1995, MONTH=1,
                 released=10, labrador=5, scotian=1),
            dict(DATE=datetime(1995, 3, 1), YEAR=1995, MONTH=3,
                 released=10, labrador=5, scotian=1),
        ])
        out = track.to_monthly(per_release)
        assert out["MONTH"].tolist() == [1, 3]

    def test_output_is_sorted_by_date(self):
        per_release = pd.DataFrame([
            dict(DATE=datetime(1996, 2, 1), YEAR=1996, MONTH=2,
                 released=10, labrador=1, scotian=1),
            dict(DATE=datetime(1995, 12, 1), YEAR=1995, MONTH=12,
                 released=10, labrador=1, scotian=1),
        ])
        out = track.to_monthly(per_release)
        assert out[["YEAR", "MONTH"]].values.tolist() == [[1995, 12], [1996, 2]]

    def test_index_sign_follows_the_definition(self):
        """More to Labrador than Scotian is positive: high retroflection."""
        retroflecting = pd.DataFrame([
            dict(DATE=datetime(1995, 1, 1), YEAR=1995, MONTH=1,
                 released=100, labrador=80, scotian=5),
        ])
        continuing = pd.DataFrame([
            dict(DATE=datetime(1995, 1, 1), YEAR=1995, MONTH=1,
                 released=100, labrador=5, scotian=80),
        ])
        assert track.to_monthly(retroflecting).loc[0, "index"] > 0
        assert track.to_monthly(continuing).loc[0, "index"] < 0
