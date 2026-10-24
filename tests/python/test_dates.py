"""Tests for the pm.dates() Python API."""

import os

import pytest

import popmix as pm

PREFIX = "tests/golden/dates/family_chr20_22"

pytestmark = pytest.mark.skipif(
    not os.path.exists(f"{PREFIX}.geno"), reason="DATES test data not available"
)


@pytest.fixture(scope="module")
def geno():
    return pm.read_genotypes(PREFIX)


def test_unknown_compat_is_rejected(geno):
    with pytest.raises(ValueError, match="dates-v753"):
        pm.dates(geno, "Simulation", "French", "Yoruba", compat="dates-v999")


def test_compat_changes_only_mesh_tail(geno):
    """compat="dates-v753" restores DATES's mesh quirks, which only reach long lags."""
    kwargs = dict(qbin=10, run_jackknife=False)
    corrected = pm.dates(geno, "Simulation", "French", "Yoruba", **kwargs)
    compat = pm.dates(geno, "Simulation", "French", "Yoruba", compat="dates-v753", **kwargs)
    differ = [
        a.distance_cm
        for a, b in zip(corrected.ld_bins, compat.ld_bins)
        if (a.n_pairs, a.correlation) != (b.n_pairs, b.correlation)
    ]
    assert differ, "DATES-0001 should change this config's long-range bins"
    assert min(differ) > 90.0
