"""admixture API, bar plot and command, on the 54-sample DATES data.

The estimates are checked against ADMIXTURE 1.3.0 in tests/golden/test_admixture.py.
"""

import os
import subprocess

import numpy as np
import pandas as pd
import popmix as pm
import pytest

SMALL = "tests/golden/dates/family_chr20_22"
POPS = ["French", "Yoruba"]

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def test_q_rows_sum_to_one_and_p_is_bounded():
    r = pm.admixture(SMALL, 2, POPS)
    q = r.q_matrix()
    np.testing.assert_allclose(q.sum(axis=1), 1.0, rtol=1e-12)
    assert q.min() >= 1e-5 - 1e-15 and r.P.min() >= 1e-5 and r.P.max() <= 1 - 1e-5
    assert len(r.snps) == r.P.shape[0]


def test_the_two_populations_are_the_two_clusters():
    r = pm.admixture(SMALL, 2, POPS)
    french = r.Q[r.Q["pop"] == "French"][["Q1", "Q2"]].mean()
    yoruba = r.Q[r.Q["pop"] == "Yoruba"][["Q1", "Q2"]].mean()
    assert abs(french["Q1"] - yoruba["Q1"]) > 0.9


def test_projection_of_a_fit_sample_gives_its_q():
    fit = pm.admixture(SMALL, 2, POPS, tolerance=1e-10)
    proj = pm.admixture(SMALL, 2, ["Yoruba"], fixed_p=fit.P, tolerance=1e-10)
    q_fit = fit.Q.set_index("id").loc[proj.Q["id"], ["Q1", "Q2"]].to_numpy()
    np.testing.assert_allclose(proj.q_matrix(), q_fit, atol=1e-6)


def test_supervised_needs_k_references():
    with pytest.raises(RuntimeError, match="reference populations"):
        pm.admixture(SMALL, 3, POPS, supervised=POPS)


def test_loglik_matches_the_fit():
    r = pm.admixture(SMALL, 2, POPS)
    assert pm.admixture_loglik(SMALL, r.q_matrix(), r.P, POPS) == pytest.approx(r.loglik, rel=1e-12)


def test_bar_plot(tmp_path):
    pytest.importorskip("matplotlib")
    from popmix import plot

    r = pm.admixture(SMALL, 2, POPS)
    plot.admixture(r, path=tmp_path / "q.png")
    assert (tmp_path / "q.png").read_bytes()[:4] == b"\x89PNG"


def test_admixture_command(tmp_path):
    out = tmp_path / "run"
    r = subprocess.run(["popmix", "admixture", "--input", SMALL, "--populations", "French,Yoruba", "-k", "2",
                        "--cv", "3", "--output-prefix", str(out)], check=True, capture_output=True, text=True)
    q = pd.read_csv(f"{out}.2.Q.tsv", sep="\t")
    api = pm.admixture(SMALL, 2, POPS)
    np.testing.assert_allclose(q[["Q1", "Q2"]], api.Q[["Q1", "Q2"]], rtol=1e-11)
    assert np.loadtxt(f"{out}.2.P").shape == api.P.shape
    assert "CV error (K=2)" in r.stdout
