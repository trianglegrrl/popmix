"""Identities of f3 and f4 when a population repeats, for any populations
(TASKS.md K4). A population's f2 with itself is 0, so f3(A; A, B) = 0,
f3(A; B, B) = f2(A, B), f4(A, A; B, C) = 0, f4(A, B; A, C) = f3(A; B, C)
and f4(A, B; A, B) = f2(A, B), block by block, so the jackknife SEs agree
too. tests/golden/test_repeated_pops_admixtools2.py checks admixtools2 gives
these numbers.
"""

from __future__ import annotations

import os

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")
SETTINGS = settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
pops = st.sampled_from(POPS)
REL = 1e-10


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(SMALL, POPS)


def same(got: dict, want: dict) -> None:
    assert got["est"] == pytest.approx(want["est"], rel=REL, abs=1e-15)
    assert got["se"] == pytest.approx(want["se"], rel=REL, abs=1e-15)


@SETTINGS
@given(a=pops, b=pops)
def test_f3_with_its_target_as_a_source_is_zero(f2, a, b):
    assert f2.f3(a, a, b)["est"] == 0


@SETTINGS
@given(a=pops, b=pops)
def test_f3_with_one_source_twice_is_f2(f2, a, b):
    same(f2.f3(a, b, b), f2.f2(a, b))


@SETTINGS
@given(a=pops, b=pops, c=pops)
def test_f4_with_a_repeated_first_pair_is_zero(f2, a, b, c):
    assert f2.f4(a, a, b, c)["est"] == 0


@SETTINGS
@given(a=pops, b=pops, c=pops)
def test_f4_sharing_its_first_population_is_f3(f2, a, b, c):
    same(f2.f4(a, b, a, c), f2.f3(a, b, c))


@SETTINGS
@given(a=pops, b=pops)
def test_f4_of_a_pair_with_itself_is_f2(f2, a, b):
    same(f2.f4(a, b, a, b), f2.f2(a, b))
