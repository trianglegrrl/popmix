"""Agreement checks of the benchmark cases (benchmarks/cases.py)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2] / "benchmarks"))

import cases  # noqa: E402


def test_f2_mismatch_is_per_entry():
    """An exact zero in the reference is not an infinite relative error."""
    want = np.array([0.0, 1.0, 2.0])
    got = np.array([1e-16, 1.0 + 5e-11, 2.0])
    assert cases._f2_mismatches(got, want) == 0


def test_f2_mismatch_counts_entries_off_both_tolerances():
    want = np.array([0.0, 1.0, 0.94664])
    got = np.array([1e-14, 1.0 + 1e-9, 0.9])
    assert cases._f2_mismatches(got, want) == 3


def test_f2_mismatch_needs_matching_nans():
    assert cases._f2_mismatches(np.array([np.nan, 1.0]), np.array([np.nan, 1.0])) == 0
    assert cases._f2_mismatches(np.array([np.nan, 1.0]), np.array([0.5, 1.0])) == 1


def test_f2_case_reads_ploidy_as_admixtools2():
    """Diploid genomes admixtools2 takes as pseudo-haploid (DEVIATIONS.md ADMIXTOOLS2-0007) would differ."""
    cmd = cases.F2_CASE.popmix(Path("in"), Path("out"), 4)
    assert cmd[cmd.index("--pseudohaploid") + 1] == "admixtools2"
