"""Golden master tests: popmix smartpca vs EIGENSOFT 7.2.1 smartpca.

The reference (tests/golden/eigensoft/v7.2.1/<config>/pca.json) is smartpca
7.2.1 built to print full precision (build_smartpca.sh), checked against the
stock binary to its printed precision (generate_pca.py). Every sample outside
the config's PCA and projected populations is labeled Ignore there; popmix is
given the same populations.

Stages:

- snps: SNPs used in the final PCA and SNPs dropped as monomorphic (exact)
- outliers: samples removed, the round and eigenvector, and their scores
- eigenvalues: every eigenvalue
- coords: every sample's coordinates (PCA and projected), each eigenvector up
  to sign (smartpca's sign is LAPACK's arbitrary choice)
- tracy_widom: TW statistic, p-value and effective number of markers, with
  compat="eigensoft-7.2.1" (EIGENSOFT-0001 changes p beyond the table)
"""

import json
import os
from pathlib import Path

import numpy as np
import popmix as pm
import pytest

REF_ROOT = Path(__file__).parent / "eigensoft" / "v7.2.1"
REPO = Path(__file__).parents[2]

EIG_RTOL = 1e-9
COORD_ATOL = 1e-9
TW_RTOL = 1e-6
P_RTOL = 1e-4

CONFIGS = ["family_all", "family_lsqproject", "family_project_no_outliers"] + [
    pytest.param(c, marks=pytest.mark.slow) for c in ["aadr_1kg_lsqproject", "aadr_1kg_project"]
]

_results: dict = {}


def _run(name: str):
    if name not in _results:
        ref = json.loads((REF_ROOT / name / "pca.json").read_text())
        prefix = ref["prefix"] if os.path.isabs(ref["prefix"]) else str(REPO / ref["prefix"])
        if not os.path.exists(prefix + ".geno"):
            pytest.skip(f"genotype data not available: {prefix}")
        p = ref["params"]
        ours = pm.smartpca(prefix, ref["pca_pops"], ref["project_pops"],
                           n_components=p.get("numoutevec", 10),
                           outlier_iterations=p.get("numoutlieriter", 5),
                           lsqproject=p.get("lsqproject", "NO") == "YES",
                           compat="eigensoft-7.2.1")
        _results[name] = (ref, ours)
    return _results[name]


@pytest.mark.parametrize("name", CONFIGS)
def test_snps(name):
    ref, ours = _run(name)
    assert ours.n_snps_used == ref["n_snps_used"]
    assert ours.n_snps_killed == ref["n_snps_killed"]
    assert ours.n_pca_samples == ref["n_samples_used"]


@pytest.mark.parametrize("name", CONFIGS)
def test_outliers(name):
    ref, ours = _run(name)
    theirs = ref["outliers"]
    assert ours.outliers["id"].tolist() == [o["id"] for o in theirs]
    assert ours.outliers["iteration"].tolist() == [o["iter"] for o in theirs]
    # smartpca prints the eigenvector from 0
    assert (ours.outliers["pc"] - 1).tolist() == [o["evec"] for o in theirs]
    np.testing.assert_allclose(ours.outliers["score"].abs().to_numpy(dtype=float),
                               np.array([abs(o["sigmage"]) for o in theirs], dtype=float), rtol=1e-6)


@pytest.mark.parametrize("name", CONFIGS)
def test_eigenvalues(name):
    ref, ours = _run(name)
    theirs = np.asarray(ref["eigenvalues"])
    assert len(ours.eigenvalues) == len(theirs)
    np.testing.assert_allclose(ours.eigenvalues, theirs, rtol=EIG_RTOL, atol=1e-12)


@pytest.mark.parametrize("name", CONFIGS)
def test_coords(name):
    ref, ours = _run(name)
    theirs = ref["samples"]
    assert ours.coords["id"].tolist() == [s["id"] for s in theirs]
    assert ours.coords["pop"].tolist() == [s["pop"] for s in theirs]
    t = np.array([s["coords"] for s in theirs])
    o = ours.coords[[f"PC{j + 1}" for j in range(t.shape[1])]].to_numpy()
    for j in range(t.shape[1]):
        sign = np.sign(np.dot(o[:, j], t[:, j]))
        np.testing.assert_allclose(sign * o[:, j], t[:, j], rtol=0, atol=COORD_ATOL,
                                   err_msg=f"PC{j + 1}")


@pytest.mark.parametrize("name", CONFIGS)
def test_tracy_widom(name):
    ref, ours = _run(name)
    theirs = ref["tracy_widom"]
    assert len(ours.tracy_widom) == len(theirs)
    for col, rtol in [("eigenvalue", EIG_RTOL), ("twstat", TW_RTOL), ("p", P_RTOL), ("neff", TW_RTOL)]:
        want = np.array([np.nan if r[col if col != "eigenvalue" else "eigenvalue"] is None else r[col]
                         for r in theirs], dtype=float)
        np.testing.assert_allclose(ours.tracy_widom[col], want, rtol=rtol, atol=0, err_msg=col)


def test_corrected_p_values_differ_only_beyond_the_table():
    ref, _ = _run("family_all")
    prefix = str(REPO / ref["prefix"])
    corrected = pm.smartpca(prefix, ref["pca_pops"], n_components=4)
    beyond = np.array([r["twstat"] is not None and r["twstat"] > 8 for r in ref["tracy_widom"]])
    assert beyond.any()
    want = np.array([np.nan if r["p"] is None else r["p"] for r in ref["tracy_widom"]], dtype=float)
    got = corrected.tracy_widom["p"].to_numpy()
    np.testing.assert_allclose(got[~beyond], want[~beyond], rtol=P_RTOL)
    tw = corrected.tracy_widom["twstat"].to_numpy()[beyond]
    # the tail is the density over sqrt(tw)
    np.testing.assert_allclose(got[beyond], want[beyond] / np.sqrt(tw), rtol=1e-10)
