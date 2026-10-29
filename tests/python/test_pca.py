"""smartpca API, plots and the pca command, on the 50-sample family data.

The numbers are checked against EIGENSOFT smartpca in
tests/golden/test_pca_eigensoft.py.
"""

import os
import subprocess

import numpy as np
import pandas as pd
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def test_default_pca_uses_every_population():
    r = pm.smartpca(SMALL, n_components=3)
    assert r.n_pca_samples == 50
    assert set(r.coords["role"]) == {"pca"}
    assert list(r.coords.columns) == ["id", "pop", "role", "PC1", "PC2", "PC3"]
    assert r.eigenvalues.sum() == pytest.approx(49.0)


def test_projected_samples_are_marked():
    r = pm.smartpca(SMALL, ["French", "Yoruba"], ["Simulation.1"], n_components=2)
    roles = r.coords.set_index("id")["role"]
    assert (roles[r.coords["pop"].values == "Simulation.1"] == "project").all()
    assert r.n_pca_samples == 45


def test_a_population_cannot_be_both_pca_and_projected():
    with pytest.raises(ValueError, match="both a PCA and a projected population"):
        pm.smartpca(SMALL, ["French", "Yoruba"], ["Yoruba"])


def test_unknown_compat_is_an_error():
    with pytest.raises(ValueError, match="unknown compat"):
        pm.smartpca(SMALL, compat="eigensoft-6")


def test_pca_plot_writes_an_image(tmp_path):
    pytest.importorskip("matplotlib")
    from popmix import plot

    r = pm.smartpca(SMALL, ["French", "Yoruba"], ["Simulation", "Simulation.1"], n_components=2)
    out = tmp_path / "pca.png"
    plot.pca(r, path=out)
    assert out.read_bytes()[:4] == b"\x89PNG"


def test_pca_command_writes_tables_matching_the_api(tmp_path):
    out = tmp_path / "run"
    subprocess.run(["popmix", "pca", "--input", SMALL, "--pca-pops", "French,Yoruba",
                    "--project-pops", "Simulation,Simulation.1", "--components", "3", "--lsqproject",
                    "--output-prefix", str(out)], check=True, capture_output=True, text=True)
    api = pm.smartpca(SMALL, ["French", "Yoruba"], ["Simulation", "Simulation.1"], n_components=3,
                      lsqproject=True)
    coords = pd.read_csv(f"{out}.coords.tsv", sep="\t")
    assert coords["id"].tolist() == api.coords["id"].tolist()
    np.testing.assert_allclose(coords[["PC1", "PC2", "PC3"]], api.coords[["PC1", "PC2", "PC3"]], rtol=1e-11,
                               atol=1e-14)
    evals = np.loadtxt(f"{out}.eval")
    # the last eigenvalue is 0 up to rounding (centering removes one dimension)
    np.testing.assert_allclose(evals, api.eigenvalues, rtol=1e-11, atol=1e-12)
    tw = pd.read_csv(f"{out}.tracy_widom.tsv", sep="\t")
    assert len(tw) == len(api.tracy_widom)
    assert pd.read_csv(f"{out}.outliers.tsv", sep="\t").columns.tolist() == ["id", "pop", "iteration", "pc",
                                                                             "score"]


def test_pca_command_plot(tmp_path):
    pytest.importorskip("matplotlib")
    out = tmp_path / "run"
    subprocess.run(["popmix", "pca", "--input", SMALL, "--components", "2", "--output-prefix", str(out),
                    "--plot", str(tmp_path / "pca.png")], check=True, capture_output=True, text=True)
    assert (tmp_path / "pca.png").read_bytes()[:4] == b"\x89PNG"
