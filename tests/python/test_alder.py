"""alder API, plot and command, on the 54-sample DATES data.

The curves and fits are checked against ALDER 1.03 in tests/golden/test_alder.py.
"""

import os
import subprocess

import numpy as np
import pandas as pd
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def test_two_references_give_the_test_and_both_one_reference_curves():
    r = pm.alder(SMALL, "Simulation.1", ["French", "Yoruba"])
    assert set(r.oneref_curves) == {"French", "Yoruba"}
    assert r.test is not None and r.mixture_fraction is None
    assert len(r.curve.fits) == 5
    assert r.curve.fit["mindis_cm"] == pytest.approx(100 * r.curve.fit_start)
    assert np.isinf(r.curve.curve["d_cm"].iloc[-1])
    assert r.chromosomes == ["20", "21", "22"]


def test_one_reference_gives_a_mixture_fraction_bound():
    r = pm.alder(SMALL, "Simulation.1", ["French"])
    assert r.test is None and 0 < r.mixture_fraction[0] < 1


def test_mindis_skips_the_ld_correlation():
    r = pm.alder(SMALL, "Simulation.1", ["French", "Yoruba"], mindis=0.01)
    assert r.ld_corr_stops == {"French": None, "Yoruba": None}
    assert r.curve.fit_start == 0.01


def test_errors():
    with pytest.raises(ValueError, match="cannot also be a reference"):
        pm.alder(SMALL, "French", ["French", "Yoruba"])
    with pytest.raises(RuntimeError, match="one or two reference"):
        pm.alder(SMALL, "Simulation.1", ["French", "Yoruba", "Simulation"])


def test_plot(tmp_path):
    pytest.importorskip("matplotlib")
    from popmix import plot

    plot.alder_curve(pm.alder(SMALL, "Simulation.1", ["French", "Yoruba"]), path=tmp_path / "c.png")
    assert (tmp_path / "c.png").read_bytes()[:4] == b"\x89PNG"


def test_alder_command(tmp_path):
    out = tmp_path / "run"
    r = subprocess.run(["popmix", "alder", "--input", SMALL, "--test", "Simulation.1", "--refs", "French,Yoruba",
                        "--output-prefix", str(out)], check=True, capture_output=True, text=True)
    assert "admixture test FAILS" in r.stdout
    fits = pd.read_csv(f"{out}.fits.tsv", sep="\t")
    api = pm.alder(SMALL, "Simulation.1", ["French", "Yoruba"])
    np.testing.assert_allclose(fits["decay"], api.curve.fits["decay"], rtol=1e-11)
