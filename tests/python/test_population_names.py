"""A population given as a string is one population, as in admixtools2,
everywhere a function takes populations: ``list("CEU.DG")`` would make it
six one-character names."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

import popmix as pm
from popmix.kinship import king_tsv

SMALL = "tests/golden/dates/family_chr20_22"
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]
RIGHT = ["Yoruba", "Simulation", "Simulation.1"]

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(SMALL, POPS)


def test_qpadm_with_one_source(f2):
    assert pm.qpadm(f2, "French", RIGHT[:2], "Simulation.1").chisq == \
        pm.qpadm(f2, ["French"], RIGHT[:2], "Simulation.1").chisq


def test_qpadm_on_genotypes_with_one_source():
    assert pm.qpadm(SMALL, "French", RIGHT[:2], "Simulation.1").chisq == \
        pm.qpadm(SMALL, ["French"], RIGHT[:2], "Simulation.1").chisq


def test_king():
    pd.testing.assert_frame_equal(pm.king(SMALL, "French"), pm.king(SMALL, ["French"]))


def test_king_tsv():
    assert king_tsv(SMALL, "French") == king_tsv(SMALL, ["French"])


def test_read_kinship():
    pd.testing.assert_frame_equal(pm.read_kinship(SMALL, "French"), pm.read_kinship(SMALL, ["French"]))


def test_sample_coverage():
    pd.testing.assert_frame_equal(pm.sample_coverage(SMALL, "Yoruba"), pm.sample_coverage(SMALL, ["Yoruba"]))


def test_sfs():
    np.testing.assert_array_equal(pm.sfs(SMALL, "Yoruba").data, pm.sfs(SMALL, ["Yoruba"]).data)


def test_joint_sfs():
    pd.testing.assert_frame_equal(pm.joint_sfs(SMALL, "Yoruba"), pm.joint_sfs(SMALL, ["Yoruba"]))


def test_outgroup_f3_matrix(f2):
    pd.testing.assert_frame_equal(pm.outgroup_f3_matrix(f2, "Yoruba", "French"),
                                  pm.outgroup_f3_matrix(f2, "Yoruba", ["French"]))


def test_the_bindings_refuse_a_string():
    # PyO3 will not take a str as a list, so the Rust functions never split one
    with pytest.raises(TypeError, match="Can't extract `str` to `Vec`"):
        pm.extract_f2(SMALL, "French")


def test_rolloff_names_the_population_it_lacks():
    with pytest.raises(Exception, match="two") as e:
        pm.rolloff(SMALL, "Simulation.1", "French")
    assert "'F'" not in str(e.value)
