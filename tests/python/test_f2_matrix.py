"""f2_matrix() holds the f2() estimates (length-weighted block jackknife)."""

import os

import numpy as np
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def test_f2_matrix_equals_f2():
    f2 = pm.extract_f2(SMALL, POPS)
    assert len(set(f2.block_lengths)) > 1
    pops = list(f2.populations)
    matrix = f2.f2_matrix()
    expected = np.array([[0.0 if a == b else f2.f2(a, b)["estimate"] for b in pops] for a in pops])
    np.testing.assert_array_equal(matrix, expected)
