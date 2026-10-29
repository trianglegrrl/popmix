"""rolloff API and command, on the 50-sample family data.

The bins are checked against ROLLOFF (ADMIXTOOLS 7.0.2) in tests/golden/test_rolloff.py.
"""

import os
import subprocess

import numpy as np
import pandas as pd
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def test_bins_are_labeled_by_their_upper_edge():
    r = pm.rolloff(SMALL, "Simulation.1", ["French", "Yoruba"], binsize=0.001, maxdis=0.02)
    assert r.bins["d_cm"].iloc[0] == pytest.approx(0.1)
    assert len(r.bins) == 20
    assert len(r.jackknife_bins) == len(r.chromosomes) == 3


def test_a_chromosome_subset(tmp_path):
    r = pm.rolloff(SMALL, "Simulation.1", ["French", "Yoruba"], maxdis=0.02, chromosomes=[20, 21])
    assert r.chromosomes == ["20", "21"]


def test_errors():
    with pytest.raises(ValueError, match="two reference"):
        pm.rolloff(SMALL, "Simulation.1", ["French"])
    with pytest.raises(ValueError, match="cannot also be a reference"):
        pm.rolloff(SMALL, "French", ["French", "Yoruba"])


def test_rolloff_command(tmp_path):
    out = tmp_path / "run"
    subprocess.run(["popmix", "rolloff", "--input", SMALL, "--admixed", "Simulation.1", "--sources", "French,Yoruba",
                    "--maxdis", "0.02", "--output-prefix", str(out)], check=True, capture_output=True, text=True)
    bins = pd.read_csv(f"{out}.bins.tsv", sep="\t")
    api = pm.rolloff(SMALL, "Simulation.1", ["French", "Yoruba"], maxdis=0.02)
    np.testing.assert_allclose(bins["corr"], api.bins["corr"], rtol=1e-11, atol=1e-15)
