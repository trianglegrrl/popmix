"""Tests for the qpadm() compat argument (DEVIATIONS.md, ADMIXTOOLS2-0001)."""

import os

import pytest

import popmix as pm

PREFIX = "tests/golden/dates/family_chr20_22"

pytestmark = pytest.mark.skipif(
    not os.path.exists(f"{PREFIX}.geno"), reason="test data not available"
)


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(PREFIX, populations=["Simulation", "French", "Yoruba", "Simulation.1"],
                         auto_only=False)


def test_unknown_compat_is_rejected(f2):
    with pytest.raises(ValueError, match="admixtools2-2.0.10"):
        pm.qpadm(f2, ["French"], ["Yoruba", "Simulation.1"], "Simulation", compat="admixtools2-9")


def test_compat_leaves_a_one_source_model_unchanged(f2):
    """With one source there is no nested path, so only `best` on the full
    model can differ (popmix puts it on the path; admixtools2 marks nothing)."""
    args = (f2, ["French"], ["Yoruba", "Simulation.1"], "Simulation")
    fixed, compat = pm.qpadm(*args), pm.qpadm(*args, compat="admixtools2-2.0.10")
    assert fixed.weights == compat.weights
    assert fixed.chisq == compat.chisq
    assert [(x.pattern, x.best) for x in fixed.popdrop] == [("0", True)]
    assert [(x.pattern, x.best) for x in compat.popdrop] == [("0", False)]
