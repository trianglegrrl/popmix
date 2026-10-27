"""Golden master tests: f3 and f4 from f2 blocks with a repeated population.

family_4pop/repeated_pops.json (generate_repeated_pops.R) holds admixtools2
f3() and f4() on the small family data's f2 blocks for combinations where a
population appears twice. admixtools2 sets a population's f2 with itself to
0, so f3(A; A, B) = 0, f4(A, B; A, C) = f3(A; B, C), and so on.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

HERE = Path(__file__).parent
REF = HERE / "admixtools2" / "v2.0.10" / "family_4pop" / "repeated_pops.json"
SMALL = HERE / "dates" / "family_chr20_22"
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(str(SMALL), POPS)


def _cases(stat: str):
    ref = json.loads(REF.read_text())
    return [pytest.param(r, id=",".join(r[k] for k in r if k.startswith("pop"))) for r in ref[stat]]


def _check(got: dict, want: dict) -> None:
    np.testing.assert_allclose(got["est"], float(want["est"]), rtol=1e-10, atol=1e-15)
    np.testing.assert_allclose(got["se"], float(want["se"]), rtol=1e-6, atol=1e-15)


@pytest.mark.parametrize("want", _cases("f3"))
def test_f3_with_a_repeated_population(f2, want):
    _check(f2.f3(want["pop1"], want["pop2"], want["pop3"]), want)


@pytest.mark.parametrize("want", _cases("f4"))
def test_f4_with_a_repeated_population(f2, want):
    _check(f2.f4(want["pop1"], want["pop2"], want["pop3"], want["pop4"]), want)
