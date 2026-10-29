"""Golden master tests: popmix admixture vs ADMIXTURE 1.3.0.

The reference (tests/golden/admixture/v1.3.0/<config>/admixture.json, from
generate_admixture.py) is ADMIXTURE run with --seed=1 and -C=1e-10 on PLINK
files written by EIGENSOFT convertf. ADMIXTURE prints Q and P to 6 decimals
and the log-likelihood to 6 decimals.

- loglik: popmix's log-likelihood at ADMIXTURE's Q and P is ADMIXTURE's, to
  the precision lost by printing them (the printed Q rows renormalized, as
  each sums to 1 in ADMIXTURE only before rounding): within 4 SDs of the
  change printing Q and P to 6 decimals makes (printing_sd; 1.6e-3 at K = 3
  on the family data, whose rare alleles put many P near 0)
- unsupervised: popmix's maximum is at least ADMIXTURE's, and where the
  maximum is clearly identified (K = 2 on the small data, the three 1000
  Genomes populations at K = 3) Q agrees up to the order of the clusters
  (to 1e-4, and to 5e-3 on AADR, where ADMIXTURE stops 0.49 log-likelihood
  units short). At
  K = 3 the small data (two populations and admixed samples) has no clear
  third component: ADMIXTURE stops at a local maximum about 160
  log-likelihood units below popmix's, with a different Q
- projection: Q with ADMIXTURE's P fixed (a convex problem, one maximum)
- supervised: popmix's maximum is at least ADMIXTURE's. On the family data
  ADMIXTURE stops 1.15 log-likelihood units below it, so the two Q differ
  (by up to 3.4e-4); its Q is checked instead as the maximum given its own
  P, which projection (convex in Q) finds: the admixed samples' Q to
  ADMIXTURE's printed precision
- cross-validation: popmix's folds are its own, so the error agrees only to
  the variation between random folds, and only where the maximum is
  identified: at K = 3 on the small data each seed's folds end at another
  local maximum (popmix's CV error 0.564 to 0.594 over seeds 1-3)
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


PRINTED = 5e-7  # half a unit of the 6th decimal


def printing_sd(prefix: str, ref: dict, q: np.ndarray, p: np.ndarray) -> float:
    """SD of the change in the log-likelihood from printing Q and P.

    Each printed value is within PRINTED of ADMIXTURE's, an error uniform on
    (-PRINTED, PRINTED) (variance PRINTED^2 / 3), except values at the
    bounds 1e-5 and 1 - 1e-5, which are exact. To first order the
    log-likelihood moves by the gradient times those errors: the SD is
    sqrt(PRINTED^2 / 3 * sum of the squared gradient over the values not at
    a bound). The gradient of sum_i h log f + (2 - h) log(1 - f), f = Q P',
    in P_jk is sum_i q_ik r_ij and in Q_ik sum_j p_jk r_ij, r = h / f -
    (2 - h) / (1 - f).
    """
    geno = pm.read_genotypes(prefix)
    col = {smp.id: n for n, smp in enumerate(geno.samples)}
    data = geno.to_numpy()[:, [col[i] for i in ref["ids"]]]
    chroms = _chroms(ref)
    wanted = None if chroms is None else {str(c) for c in chroms}
    keep = np.array([wanted is None or str(snp.chromosome) in wanted for snp in geno.snps])
    data = data[keep]
    data = data[(data < 3).any(axis=1)]  # SNPs with no call are left out
    assert data.shape[0] == p.shape[0]
    called = data < 3
    h = np.where(called, 2.0 - data, 0.0)  # copies of the allele P is for
    f = p @ q.T
    r = np.where(called, h / f - (2.0 - h) / (1.0 - f), 0.0)
    printed_p, printed_q = np.array(ref["P"]), np.array(ref["Q"])
    free_p = (printed_p > 1e-5) & (printed_p < 1 - 1e-5)
    free_q = (printed_q > 1e-5) & (printed_q < 1 - 1e-5)
    squares = ((r @ q)[free_p] ** 2).sum() + ((r.T @ p)[free_q] ** 2).sum()
    return float(np.sqrt(PRINTED ** 2 / 3 * squares))


@pytest.mark.parametrize("name", UNSUPERVISED + ["family_supervised"])
def test_loglik_at_admixture_estimates(name):
    ref, prefix = _ref(name)
    q = np.array(ref["Q"])
    q /= q.sum(axis=1, keepdims=True)
    p = np.array(ref["P"])
    ll = pm.admixture_loglik(prefix, q, p, ref["pops"], chromosomes=_chroms(ref))
    # 4 SDs of the printing error, plus the printed log-likelihood's rounding
    assert abs(ll - ref["loglik"]) <= 4 * printing_sd(prefix, ref, q, p) + PRINTED


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
    assert r.Q["id"].tolist() == ref["ids"]
    assert r.loglik >= ref["loglik"] - 1e-6
    # ADMIXTURE's Q is the maximum at its P: projection on that P finds it
    at_p = pm.admixture(prefix, ref["k"], ref["pops"], fixed_p=np.array(ref["P"]), tolerance=TOL)
    admixed = (~r.Q["pop"].isin(ref["references"])).to_numpy()
    np.testing.assert_allclose(at_p.q_matrix()[admixed], np.array(ref["Q"])[admixed], atol=2e-6)


@pytest.mark.parametrize("name", [n for n in UNSUPERVISED if n != "family_k3"])
def test_cross_validation(name):
    ref, prefix = _ref(name)
    r = pm.admixture(prefix, ref["k"], ref["pops"], chromosomes=_chroms(ref), cv=5, tolerance=1e-6)
    print(f"{name}: CV error {r.cv_error:.5f}, ADMIXTURE {ref['cv_error']:.5f}")
    assert r.cv_error == pytest.approx(ref["cv_error"], rel=0.02)
