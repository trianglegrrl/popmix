"""Bin options that cannot work are a ValueError, never a panic or an abort.

rolloff(binsize=0) used to overflow an allocation (PanicException) and a tiny
binsize to abort the interpreter allocating its bins.
"""

import math
import os

import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")

BAD_BINSIZES = [0.0, -0.001, math.nan, math.inf, 1e-12]


@pytest.mark.parametrize("binsize", BAD_BINSIZES)
def test_rolloff_binsize(binsize):
    with pytest.raises(ValueError, match="binsize"):
        pm.rolloff(SMALL, "Simulation.1", ["French", "Yoruba"], binsize=binsize, maxdis=0.02)


@pytest.mark.parametrize("maxdis", [0.0, math.nan, -1.0])
def test_rolloff_maxdis(maxdis):
    with pytest.raises(ValueError, match="maxdis"):
        pm.rolloff(SMALL, "Simulation.1", ["French", "Yoruba"], maxdis=maxdis)


@pytest.mark.parametrize("binsize", BAD_BINSIZES)
def test_alder_binsize(binsize):
    with pytest.raises(ValueError, match="binsize"):
        pm.alder(SMALL, "Simulation.1", ["French", "Yoruba"], binsize=binsize)


@pytest.fixture(scope="module")
def geno():
    return pm.read_genotypes(SMALL)


@pytest.mark.parametrize("width", [0.0, -0.1, math.nan, 1e-12])
def test_dates_bin_width(geno, width):
    with pytest.raises(ValueError, match="bin_width_cm"):
        pm.dates(geno, "Simulation.1", "French", "Yoruba", bin_width_cm=width, run_jackknife=False)


def test_dates_qbin(geno):
    with pytest.raises(ValueError, match="qbin"):
        pm.dates(geno, "Simulation.1", "French", "Yoruba", qbin=10**9, run_jackknife=False)


@pytest.mark.parametrize("years", [0.0, math.nan])
def test_dates_generation_time(geno, years):
    with pytest.raises(ValueError, match="generation_time"):
        pm.dates(geno, "Simulation.1", "French", "Yoruba", generation_time=years, run_jackknife=False)
