"""Golden master tests: popmix's qpAdm model-set functions vs admixtools2.

The reference (tests/golden/admixtools2/v2.0.10/dg_11pop/qpadm_models.json,
from generate_qpadm_models.R) holds, on the dg_11pop f2 blocks:

- qpadm_p for full-rank models (with and without constraints), a rank-0
  model and one without a target, and test_cladality
- qpwave_pairs for five left populations
- qpadm_rotate and rotate_models: five populations rotated between left and
  right with two fixed right populations
- qpadm_multi: three models, summary rows and full results
- qpadm_models: the model sets for 4 and 5 populations

Statistics compare at CLAUDE.md's tolerances: chi-squared (through the
inverse f4 covariance, as standard errors; as in test_fstats_admixtools2)
rtol 1e-6, p-values rtol 1e-4, weights rtol 1e-6; ranks, degrees of freedom,
feasibility and the models themselves exactly.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

REF_DIR = Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_11pop"
REF_FILE = REF_DIR / "qpadm_models.json"
CHISQ_RTOL = 1e-6
P_RTOL = 1e-4
WEIGHT_RTOL = 1e-6

pytestmark = pytest.mark.skipif(not REF_FILE.exists(), reason="qpadm_models.json not generated")


def _ref() -> dict:
    return json.loads(REF_FILE.read_text())


def _list(x) -> list:
    return x if isinstance(x, list) else [x]


def _num(x):
    return np.nan if x in (None, "NA") else float(x)


@pytest.fixture(scope="module")
def f2_blocks():
    f2_ref = json.loads((REF_DIR / "f2.json").read_text())
    if not os.path.exists(f"{f2_ref['prefix']}.geno"):
        pytest.skip(f"genotype data not available: {f2_ref['prefix']}")
    assert f2_ref["pops"] == _ref()["pops"]
    p = f2_ref["params"]
    # admixtools2's ploidy rule (ADMIXTOOLS2-0007): Druze.DG holds a diploid
    # genome it takes as pseudo-haploid (test_fstats_admixtools2.py)
    return pm.extract_f2(f2_ref["prefix"], f2_ref["pops"], auto_only=p["auto_only"], block_size=p["blgsize"],
                         maxmiss=p["maxmiss"], poly_only=p["poly_only"] == "f2", pseudohaploid="admixtools2")


def _check_stats(ours: dict, theirs: dict, where: str):
    assert ours["f4rank"] == int(theirs["f4rank"]), where
    assert ours["dof"] == int(theirs["dof"]), where
    assert ours["chisq"] == pytest.approx(float(theirs["chisq"]), rel=CHISQ_RTOL), where
    assert _num(ours["p"]) == pytest.approx(_num(theirs["p"]), rel=P_RTOL, abs=1e-300, nan_ok=True), where
    assert ours["feasible"] == theirs["feasible"], where


def test_qpadm_p(f2_blocks):
    for case in _ref()["qpadm_p"]:
        ours = pm.qpadm_p(f2_blocks, case["left"], case["right"], case.get("target"), rnk=case.get("rnk"),
                          constrained=case.get("constrained", False), weights=True)
        where = f"{case['target']} ~ {case['left']}"
        _check_stats(ours, case, where)
        theirs = case["weights"]
        if theirs is None:
            assert ours["weights"] is None, where
        else:
            assert list(ours["weights"]) == list(theirs), where
            np.testing.assert_allclose(list(ours["weights"].values()), [float(v) for v in theirs.values()],
                                       rtol=WEIGHT_RTOL, err_msg=where)


def test_test_cladality(f2_blocks):
    t = _ref()["test_cladality"]
    ours = pm.test_cladality(f2_blocks, t["left"], t["right"])
    assert ours["dof"] == t["dof"]
    assert ours["chisq"] == pytest.approx(float(t["chisq"]), rel=CHISQ_RTOL)
    assert ours["p"] == pytest.approx(float(t["p"]), rel=P_RTOL)


def test_qpwave_pairs(f2_blocks):
    t = _ref()["qpwave_pairs"]
    ours = pm.qpwave_pairs(f2_blocks, t["left"], t["right"]).to_dict("records")
    assert [(r["pop1"], r["pop2"]) for r in ours] == [(r["pop1"], r["pop2"]) for r in t["result"]]
    for o, r in zip(ours, t["result"]):  # noqa: B905
        assert o["chisq"] == pytest.approx(float(r["chisq"]), rel=CHISQ_RTOL), (o["pop1"], o["pop2"])
        assert o["p"] == pytest.approx(float(r["p"]), rel=P_RTOL, abs=1e-300), (o["pop1"], o["pop2"])


def test_rotate_models():
    t = _ref()["rotate_models"]
    ours = pm.rotate_models(t["leftright"], t["target"], t["rightfix"]).to_dict("records")
    assert [(m["left"], m["right"], m["target"]) for m in ours] == \
        [(_list(m["left"]), _list(m["right"]), m["target"]) for m in t["models"]]


def test_qpadm_rotate(f2_blocks):
    t = _ref()["qpadm_rotate"]
    ours = pm.qpadm_rotate(f2_blocks, t["leftright"], t["target"], t["rightfix"]).to_dict("records")
    assert [(r["left"], r["right"]) for r in ours] == [(_list(r["left"]), _list(r["right"])) for r in t["result"]]
    for o, r in zip(ours, t["result"]):  # noqa: B905
        _check_stats(o, r, f"{o['left']} | {o['right']}")


def test_qpadm_multi_summary(f2_blocks):
    t = _ref()["qpadm_multi"]
    models = [{"left": _list(m["left"]), "right": _list(m["right"]), "target": m["target"]} for m in t["models"]]
    ours = pm.qpadm_multi(f2_blocks, models, full_results=False).to_dict("records")
    assert len(ours) == len(t["result"])
    for o, r in zip(ours, t["result"]):  # noqa: B905
        assert (o["left"], o["right"], o["target"]) == (_list(r["left"]), _list(r["right"]), r["target"])
        _check_stats(o, r, f"{o['target']} ~ {o['left']}")


def test_qpadm_multi_full_results(f2_blocks):
    t = _ref()["qpadm_multi"]
    models = [{"left": _list(m["left"]), "right": _list(m["right"]), "target": m["target"]} for m in t["models"]]
    ours = pm.qpadm_multi(f2_blocks, models)
    for res, full, m in zip(ours, t["full"], models):  # noqa: B905
        assert res.sources == tuple(m["left"])
        np.testing.assert_allclose(res.weights, [float(w["weight"]) for w in full["weights"]], rtol=WEIGHT_RTOL)
        # a lone source has weight 1 and SE 0, up to rounding
        np.testing.assert_allclose(res.weight_se, [float(w["se"]) for w in full["weights"]], rtol=1e-6,
                                   atol=1e-12)
        assert res.chisq == pytest.approx(float(full["popdrop"][0]["chisq"]), rel=CHISQ_RTOL)
        assert [r.pattern for r in res.popdrop] == [r["pat"] for r in full["popdrop"]]


def test_qpadm_models():
    for t in _ref()["qpadm_models"]:
        ours = pm.qpadm_models(t["pops"], allpops=t["allpops"], more_right=t["more_right"]).to_dict("records")
        assert [(m["left"], m["right"]) for m in ours] == \
            [(_list(m["left"]), _list(m["right"])) for m in t["models"]], (t["allpops"], t["more_right"])
