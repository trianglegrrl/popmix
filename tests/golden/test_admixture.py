"""Golden master tests: popmix admixture vs ADMIXTURE 1.3.0.

The reference (tests/golden/admixture/v1.3.0/<config>/admixture.json, from
generate_admixture.py) is ADMIXTURE run with --seed=1 and -C=1e-10 on PLINK
files written by EIGENSOFT convertf. ADMIXTURE prints Q and P to 6 decimals
and the log-likelihood to 6 decimals.

- loglik: popmix's log-likelihood at ADMIXTURE's Q and P is ADMIXTURE's (to
  the precision lost by printing them; the printed Q rows are renormalized,
  as each sums to 1 in ADMIXTURE only before rounding)
- unsupervised: popmix's maximum is at least ADMIXTURE's, and where the
  maximum is clearly identified (K = 2 on the small data, the three 1000
  Genomes populations at K = 3) Q agrees up to the order of the clusters
  (to 1e-4, and to 5e-3 on AADR, where ADMIXTURE stops 0.49 log-likelihood
  units short). At
  K = 3 the small data (two populations and admixed samples) has no clear
  third component: ADMIXTURE stops at a local maximum about 160
  log-likelihood units below popmix's, with a different Q
- projection: Q with ADMIXTURE's P fixed (a convex problem, one maximum)
- supervised: as unsupervised, with the reference samples fixed
- cross-validation: popmix's folds are its own, so the error agrees only to
  the variation between random folds
"""

import json
import os
from itertools import permutations
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

from .reference_paths import local, missing

REF_ROOT = Path(__file__).parent / "admixture" / "v1.3.0"
REPO = Path(__file__).parents[2]
SLOW = pytest.mark.slow
UNSUPERVISED = ["family_k2", "family_k3", pytest.param("aadr_1kg_k3", marks=SLOW)]
# Q tolerance where the maximum is identified. On aadr_1kg_k3 ADMIXTURE stops
# 0.49 log-likelihood units below popmix's maximum (out of 1.5e7), along a
# flat direction, and Q differs by up to 3e-3.
IDENTIFIED = {"family_k2": 1e-4, "aadr_1kg_k3": 5e-3}
PROJECTION = ["family_projection", pytest.param("aadr_1kg_projection", marks=SLOW)]
TOL = 1e-10


def _ref(name: str) -> tuple[dict, str]:
    ref = json.loads((REF_ROOT / name / "admixture.json").read_text(encoding="utf-8"))
    prefix = local(ref["prefix"]) if os.path.isabs(local(ref["prefix"])) else str(REPO / local(ref["prefix"]))
    if missing(prefix + ".geno"):
        pytest.skip(f"genotype data not available: {prefix}")
    return ref, prefix


def _chroms(ref: dict):
    c = ref.get("chromosomes")
    if c is None:
        return None
    lo, hi = c.split("-")
    return list(range(int(lo), int(hi) + 1))


def _match(ours: np.ndarray, theirs: np.ndarray) -> np.ndarray:
    """ours with its columns permuted to best match theirs."""
    best = min(permutations(range(theirs.shape[1])), key=lambda p: np.abs(ours[:, p] - theirs).max())
    return ours[:, best]


@pytest.mark.parametrize("name", UNSUPERVISED + ["family_supervised"])
def test_loglik_at_admixture_estimates(name):
    ref, prefix = _ref(name)
    q = np.array(ref["Q"])
    q /= q.sum(axis=1, keepdims=True)
    ll = pm.admixture_loglik(prefix, q, np.array(ref["P"]), ref["pops"], chromosomes=_chroms(ref))
    # Q and P printed to 6 decimals move the log-likelihood by about 1e-9 of it
    assert ll == pytest.approx(ref["loglik"], rel=1e-9)


@pytest.mark.parametrize("name", UNSUPERVISED)
def test_unsupervised(name):
    ref, prefix = _ref(name)
    r = pm.admixture(prefix, ref["k"], ref["pops"], chromosomes=_chroms(ref), tolerance=TOL)
    assert r.Q["id"].tolist() == ref["ids"]
    assert r.loglik >= ref["loglik"] - 1e-6
    if name in IDENTIFIED:
        np.testing.assert_allclose(_match(r.q_matrix(), np.array(ref["Q"])), ref["Q"], atol=IDENTIFIED[name])


@pytest.mark.parametrize("name", PROJECTION)
def test_projection(name):
    ref, prefix = _ref(name)
    r = pm.admixture(prefix, ref["k"], ref["project"], fixed_p=np.array(ref["P_in"]), chromosomes=_chroms(ref),
                     snps=ref["snps"], tolerance=TOL)
    assert r.Q["id"].tolist() == ref["ids"]
    assert r.loglik >= ref["loglik"] - 1e-6
    np.testing.assert_allclose(r.q_matrix(), ref["Q"], atol=2e-6)


def test_supervised():
    ref, prefix = _ref("family_supervised")
    r = pm.admixture(prefix, ref["k"], ref["pops"], supervised=ref["references"], tolerance=TOL)
    assert r.loglik >= ref["loglik"] - 1e-6
    np.testing.assert_allclose(r.q_matrix(), ref["Q"], atol=1e-4)


@pytest.mark.parametrize("name", UNSUPERVISED)
def test_cross_validation(name):
    ref, prefix = _ref(name)
    r = pm.admixture(prefix, ref["k"], ref["pops"], chromosomes=_chroms(ref), cv=5, tolerance=1e-6)
    print(f"{name}: CV error {r.cv_error:.5f}, ADMIXTURE {ref['cv_error']:.5f}")
    assert r.cv_error == pytest.approx(ref["cv_error"], rel=0.02)
