"""DATES on the chr20-22 subset against DATES v753's LD table, bin by bin.

The full stage-by-stage comparison with DATES v753 (LD table, fit,
jackknife, the corrected default) is test_dates_v753.py.
"""

import inspect
from pathlib import Path

import pytest

import popmix as pm

from .reference_paths import missing

REF = Path(__file__).parent / "dates" / "v753" / "chr20_22_q0"


def test_dates_function_signature():
    """Test that dates function has correct signature."""
    sig = inspect.signature(pm.dates)
    params = list(sig.parameters.keys())

    assert params[:3] == ["prefix", "admixed", "sources"]
    for option in ("generation_time", "min_distance_cm", "max_distance_cm", "bin_width_cm", "jackknife"):
        assert sig.parameters[option].kind is inspect.Parameter.KEYWORD_ONLY


def test_dates_vs_reference():
    """DATES defaults (exact mode, 0.1 cM bins to 100 cM) against DATES v753's
    exact-mode LD table on the chr20-22 data (v753/chr20_22_q0/dates.out).

    Criteria:
    - Pair counts must match exactly (0 difference)
    - Correlations must match within 1e-6
    """
    data_path = "tests/golden/dates/family_chr20_22"

    if missing(f"{data_path}.geno"):
        pytest.skip("Chr20-22 subset data not available")

    ref_bins = []
    for line in (REF / "dates.out").read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            dist, _cov, _yreg, corr, n_pairs = line.split()
            ref_bins.append((float(dist), float(corr), int(n_pairs)))
    assert len(ref_bins) == 1000

    source1, source2, admixed, _ = (REF / "admixlist").read_text(encoding="utf-8").split()
    result = pm.dates(data_path, admixed, [source1, source2])
    ours = {round(b.distance_cm, 3): b for b in result.ld_bins}
    for dist_cm, corr, n_pairs in ref_bins:
        b = ours.get(round(dist_cm, 3))
        assert b is not None, f"Missing bin at {dist_cm}cM"
        assert b.n_pairs == n_pairs, f"Pair count mismatch at {dist_cm}cM: got {b.n_pairs}, expected {n_pairs}"
        assert abs(b.correlation - corr) < 1e-6, f"Correlation mismatch at {dist_cm}cM"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
