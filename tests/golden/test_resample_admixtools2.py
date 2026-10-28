"""Golden master tests: popmix resampling vs admixtools2.

The reference (tests/golden/admixtools2/v2.0.10/dg_6pop/resample.json, from
generate_resample.R) holds, on the dg_6pop f2 blocks:

- boo_list: the bootstrap block draws after set.seed(1)
- qpgraph_resample_multi for a tree and a one-admixture graph (6 bootstrap
  replicates, qpgraph seed 1), and compare_fits of their test scores
- qpgraph_resample_snps with every block left out in turn, and with 5
  bootstrap replicates
- qpadm_resample_snps (jackknife), qpwave_resample_snps (bootstrap) and
  qpdstat / qp3pop _resample_snps (jackknife)

popmix draws the same bootstrap blocks from the same seed (R's sample()), so
every replicate is compared with its counterpart: statistics at the
f-statistic tolerance, chi-squared and SEs at 1e-6, p-values at 1e-4, and
admixture weights at 1e-6. qpgraph scores and weights are optimized (from
admixtools2's starts, with R's L-BFGS-B: DEVIATIONS.md): scores at 1e-8,
weights at 1e-6 or, where the score is too flat to pin them down, within
admixtools2's spread over its starts, and out-of-sample scores (which the
fit does not minimize, so they move with the weights) at 1e-6.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

from .reference_paths import local

REF_DIR = Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_6pop"
REF_FILE = REF_DIR / "resample.json"
F_RTOL = 1e-10
SE_RTOL = 1e-6
WEIGHT_TOL = 1e-6
P_RTOL = 1e-4
# Optimized qpgraph scores: L-BFGS-B stops at a relative reduction of 2e-9
OPT_RTOL = 1e-8

pytestmark = pytest.mark.skipif(not REF_FILE.exists(), reason="resample.json not generated")


def _ref() -> dict:
    return json.loads(REF_FILE.read_text())


def _f(x) -> float:
    return np.nan if x in (None, "NA") else float(x)


def _list(x) -> list:
    return x if isinstance(x, list) else [x]


@pytest.fixture(scope="module")
def f2_blocks():
    f2_ref = json.loads((REF_DIR / "f2.json").read_text())
    if not os.path.exists(f"{local(f2_ref['prefix'])}.geno"):
        pytest.skip(f"genotype data not available: {local(f2_ref['prefix'])}")
    assert f2_ref["pops"] == _ref()["pops"]
    p = f2_ref["params"]
    return pm.extract_f2(local(f2_ref["prefix"]), f2_ref["pops"], auto_only=p["auto_only"], block_size=p["blgsize"],
                         maxmiss=p["maxmiss"], poly_only=p["poly_only"] == "f2")


def _graph(name: str) -> list:
    return [tuple(e) for e in _ref()["graphs"][name]]


def test_boo_list_draws_admixtools2s_blocks(f2_blocks):
    t = _ref()["boo_list"]
    b = pm.boo_list(f2_blocks, nboot=t["nboot"], seed=t["seed"])
    assert [list(np.add(s, 1)) for s in b["sel"]] == t["sel"]
    n = f2_blocks.n_blocks
    for boo, test, sel in zip(b["boo"], b["test"], b["sel"]):  # noqa: B905
        assert boo.n_blocks == n
        assert boo.block_lengths == [f2_blocks.block_lengths[i] for i in sel]
        assert test.n_blocks == n - len(set(sel))


def _check_qpgraph_rows(ours, theirs: list, label: str, test_score: bool = False):
    assert len(ours) == len(theirs), label
    for (_, o), r in zip(ours.iterrows(), theirs):  # noqa: B905
        where = f"{label} replicate {r['id']}"
        assert o["id"] == r["id"], where
        assert not r["error"] and o["error"] is None, where
        if r.get("sel") is not None:
            assert list(np.add(o["sel"], 1)) == _list(r["sel"]), where
        assert o["score"] == pytest.approx(float(r["score"]), rel=OPT_RTOL), where
        if test_score:
            # not minimized, so first-order in the weights' differences
            assert o["score_test"] == pytest.approx(float(r["score_test"]), rel=SE_RTOL), where
        if r["weight"] is not None:
            # where the score barely depends on it, the weight is pinned down
            # only to within admixtools2's own spread over its starts
            spread = float(r["high"]) - float(r["low"])
            assert abs(o["weights"][0] - float(r["weight"])) <= max(WEIGHT_TOL, spread), where


@pytest.mark.slow
def test_qpgraph_resample_multi(f2_blocks):
    t = _ref()["qpgraph_resample_multi"]
    names = list(t["fits"])
    fits = pm.qpgraph_resample_multi(f2_blocks, [_graph(n) for n in names], nboot=t["nboot"], boot_seed=t["seed"],
                                     seed=t["qpgraph_seed"])
    for name, ours in zip(names, fits):  # noqa: B905
        _check_qpgraph_rows(ours, t["fits"][name], name, test_score=True)


@pytest.mark.slow
def test_compare_fits_of_resampled_graphs(f2_blocks):
    t = _ref()["qpgraph_resample_multi"]
    names = list(t["fits"])
    s1, s2 = ([float(r["score_test"]) for r in t["fits"][n]] for n in names)
    ours = pm.compare_fits(s1, s2)
    for k, v in _ref()["compare_fits"].items():
        assert ours[k] == pytest.approx(float(v), rel=1e-12), k


def test_compare_fits():
    t = _ref()["compare_fits_small"]
    ours = pm.compare_fits(t["scores1"], t["scores2"])
    assert set(ours) == set(t["result"])
    for k, v in t["result"].items():
        assert ours[k] == pytest.approx(float(v), rel=1e-12), k


def test_compare_fits_needs_equal_lengths():
    with pytest.raises(ValueError, match="differ"):
        pm.compare_fits([1.0, 2.0], [1.0])


@pytest.mark.slow
def test_qpgraph_resample_snps_jackknife(f2_blocks):
    ours = pm.qpgraph_resample_snps(f2_blocks, _graph("admix"), seed=1)
    _check_qpgraph_rows(ours, _ref()["qpgraph_resample_snps_jack"], "jackknife")


def test_qpgraph_resample_snps_bootstrap(f2_blocks):
    t = _ref()["qpgraph_resample_snps_boot"]
    ours = pm.qpgraph_resample_snps(f2_blocks, _graph("admix"), boot=t["boot"], boot_seed=t["seed"], seed=1)
    _check_qpgraph_rows(ours, t["rows"], "bootstrap")


def test_qpadm_resample_snps(f2_blocks):
    t = _ref()["qpadm_resample_snps"]
    ours = pm.qpadm_resample_snps(f2_blocks, t["left"], t["right"], t["target"])
    assert len(ours) == len(t["rows"])
    for (_, o), r in zip(ours.iterrows(), t["rows"]):  # noqa: B905
        np.testing.assert_allclose(o["weights"], [float(w) for w in r["weights"]], rtol=0, atol=WEIGHT_TOL)
        assert o["chisq"] == pytest.approx(float(r["chisq"]), rel=SE_RTOL), r["id"]
        assert o["p"] == pytest.approx(float(r["p"]), rel=P_RTOL), r["id"]


def test_qpwave_resample_snps(f2_blocks):
    t = _ref()["qpwave_resample_snps"]
    ours = pm.qpwave_resample_snps(f2_blocks, t["left"], t["right"], boot=t["boot"], boot_seed=t["seed"])
    assert len(ours) == len(t["rows"])
    for (_, o), r in zip(ours.iterrows(), t["rows"]):  # noqa: B905
        assert list(np.add(o["sel"], 1)) == r["sel"]
        rd = o["rankdrop"]
        assert [x["f4rank"] for x in rd] == r["f4rank"]
        np.testing.assert_allclose([x["chisq"] for x in rd], [float(c) for c in r["chisq"]], rtol=SE_RTOL)
        np.testing.assert_allclose([x["p"] for x in rd], [float(p) for p in r["p"]], rtol=P_RTOL, atol=1e-300)


@pytest.mark.parametrize("fun", ["qpdstat_resample_snps", "qp3pop_resample_snps"])
def test_fstat_resample_snps(f2_blocks, fun):
    t = _ref()[fun]
    pops = dict(zip(["pop1", "pop2", "pop3", "pop4"], t["pops"]))  # noqa: B905
    ours = getattr(pm, fun)(f2_blocks, **pops)
    assert len(ours) == len(t["rows"])
    np.testing.assert_array_equal(ours["id"], [r["id"] for r in t["rows"]])
    np.testing.assert_allclose(ours["est"], [float(r["est"]) for r in t["rows"]], rtol=F_RTOL, atol=1e-15)
    np.testing.assert_allclose(ours["se"], [float(r["se"]) for r in t["rows"]], rtol=SE_RTOL)


INDS_FILE = REF_DIR / "resample_inds.json"


@pytest.fixture(scope="module")
def inds_ref():
    if not INDS_FILE.exists():
        pytest.skip("resample_inds.json not generated")
    ref = json.loads(INDS_FILE.read_text())
    if not os.path.exists(f"{local(ref['prefix'])}.geno"):
        pytest.skip(f"genotype data not available: {local(ref['prefix'])}")
    return ref


def _f2_args(ref: dict) -> dict:
    p = ref["params"]
    return {"auto_only": p["auto_only"], "block_size": p["blgsize"], "maxmiss": p["maxmiss"],
            "poly_only": p["poly_only"] == "f2"}


@pytest.mark.slow
def test_resample_inds(inds_ref):
    """Each statistic with each individual left out, against admixtools2 on
    f2_from_geno(inds, pops) without that individual."""
    ref = inds_ref
    args = (local(ref["prefix"]), ref["inds"], ref["pops"])
    f2a = _f2_args(ref)
    qa = pm.qpadm_resample_inds(*args, ref["qpadm"]["left"], ref["qpadm"]["right"], ref["qpadm"]["target"],
                                f2_args=f2a)
    qw = pm.qpwave_resample_inds(*args, ref["qpwave"]["left"], ref["qpwave"]["right"], f2_args=f2a)
    qd = pm.qpdstat_resample_inds(*args, *ref["qpdstat"], f2_args=f2a)
    q3 = pm.qp3pop_resample_inds(*args, *ref["qp3pop"], f2_args=f2a)
    qg = pm.qpgraph_resample_inds(*args, [tuple(e) for e in ref["graph"]], f2_args=f2a, seed=1)
    rows = ref["rows"]
    assert list(qa.ind) == [r["ind"] for r in rows] == list(qd.ind) == list(qg.ind)
    for k, r in enumerate(rows):
        where = r["ind"]
        assert qa.error[k] is None and qg.error[k] is None, where
        np.testing.assert_allclose(qa.weights[k], [float(w) for w in r["qpadm"]["weights"]], rtol=0, atol=WEIGHT_TOL)
        np.testing.assert_allclose(qa.result[k].weight_se, [float(w) for w in r["qpadm"]["se"]], rtol=SE_RTOL)
        assert qa.chisq[k] == pytest.approx(float(r["qpadm"]["chisq"]), rel=SE_RTOL), where
        rd = qw.rankdrop[k]
        assert [x["f4rank"] for x in rd] == r["qpwave"]["f4rank"], where
        np.testing.assert_allclose([x["chisq"] for x in rd], [float(c) for c in r["qpwave"]["chisq"]], rtol=SE_RTOL)
        assert qd.est[k] == pytest.approx(float(r["qpdstat"]["est"]), rel=F_RTOL), where
        assert qd.se[k] == pytest.approx(float(r["qpdstat"]["se"]), rel=SE_RTOL), where
        assert q3.est[k] == pytest.approx(float(r["qp3pop"]["est"]), rel=F_RTOL), where
        assert q3.se[k] == pytest.approx(float(r["qp3pop"]["se"]), rel=SE_RTOL), where
        assert qg.score[k] == pytest.approx(float(r["qpgraph"]["score"]), rel=OPT_RTOL), where
        spread = float(r["qpgraph"]["high"]) - float(r["qpgraph"]["low"])
        assert abs(qg.weights[k][0] - float(r["qpgraph"]["weight"])) <= max(WEIGHT_TOL, spread), where
