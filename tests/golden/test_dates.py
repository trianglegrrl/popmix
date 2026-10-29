"""DATES on the chr20-22 subset against DATES v753's LD table, bin by bin.

The full stage-by-stage comparison with DATES v753 (LD table, fit,
jackknife, the corrected default) is test_dates_v753.py.
"""

import inspect
import json

import pytest

import popmix as pm

from .reference_paths import missing


def test_dates_function_signature():
    """Test that dates function has correct signature."""
    sig = inspect.signature(pm.dates)
    params = list(sig.parameters.keys())

    assert params[:3] == ["prefix", "admixed", "sources"]
    for option in ("generation_time", "min_distance_cm", "max_distance_cm", "bin_width_cm", "jackknife"):
        assert sig.parameters[option].kind is inspect.Parameter.KEYWORD_ONLY


def test_dates_vs_reference():
    """Compare DATES output against DATES v753 reference on chr20-22 subset.

    This test validates that popmix DATES produces identical results to
    the original DATES v753 implementation. Uses the chr20-22 subset
    which runs quickly while still validating the algorithm.

    Criteria:
    - Pair counts must match exactly (0 difference)
    - Correlations must match within 1e-6
    """
    data_path = "tests/golden/dates/family_chr20_22"

    if missing(f"{data_path}.geno"):
        pytest.skip("Chr20-22 subset data not available")

    # Load reference values from DATES v753
    with open("tests/golden/dates/chr20_22_reference.json", encoding="utf-8") as f:
        ref = json.load(f)

    ref_bins = ref["ld_decay_bins"]

    # run DATES

    result = pm.dates(data_path, "Simulation", ["French", "Yoruba"])

    # Compare each bin
    total_pair_diff = 0
    max_corr_diff = 0.0

    for ref_bin in ref_bins:
        dist_cm = ref_bin["distance_cm"]
        # Find matching bin in our output
        pm_bin = None
        for b in result.ld_bins:
            if abs(b.distance_cm - dist_cm) < 0.01:
                pm_bin = b
                break

        assert pm_bin is not None, f"Missing bin at {dist_cm}cM"

        # Pair counts must match exactly
        pair_diff = pm_bin.n_pairs - ref_bin["n_pairs"]
        total_pair_diff += abs(pair_diff)

        # Correlations must match within tolerance
        corr_diff = abs(pm_bin.correlation - ref_bin["correlation"])
        max_corr_diff = max(max_corr_diff, corr_diff)

        assert pair_diff == 0, \
            f"Pair count mismatch at {dist_cm}cM: got {pm_bin.n_pairs}, expected {ref_bin['n_pairs']}"
        assert corr_diff < 1e-6, \
            f"Correlation mismatch at {dist_cm}cM: diff={corr_diff:.2e}"

    print("\nDATES vs Reference Validation:")
    print(f"  Total pair count difference: {total_pair_diff}")
    print(f"  Max correlation difference: {max_corr_diff:.2e}")
    print(f"  Coefficient α: {ref['dates_output']['coefficient_source1']}")
    print("  Status: PASS")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
