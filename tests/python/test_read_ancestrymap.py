"""An ANCESTRYMAP dataset (``.ancestrymapgeno``, ``.snp``, ``.ind``) is
read like the other formats popmix writes: found from its prefix by
read_genotypes and extract_f2."""

from __future__ import annotations

import os
import subprocess

import numpy as np
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


@pytest.fixture(scope="module")
def prefix(tmp_path_factory) -> str:
    out = str(tmp_path_factory.mktemp("am") / "data")
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", out, "-of", "ancestrymap"], check=True,
                   capture_output=True)
    return out


def test_read_genotypes(prefix):
    np.testing.assert_array_equal(pm.read_genotypes(prefix).to_numpy(), pm.read_genotypes(SMALL).to_numpy())


def test_extract_f2(prefix):
    assert pm.extract_f2(prefix, POPS).f2("French", "Yoruba") == pm.extract_f2(SMALL, POPS).f2("French", "Yoruba")
