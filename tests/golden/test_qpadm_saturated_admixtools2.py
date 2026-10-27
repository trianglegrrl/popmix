"""Golden master tests: a saturated qpAdm model is not tested.

dg_6pop/qpadm_saturated.json (generate_saturated_qpadm.R) holds admixtools2
qpadm() of French.DG = Han.DG + Yoruba.DG with right Mbuti.DG, Chimp.REF:
as many right populations as left, so the full model has 0 degrees of
freedom and admixtools2 reports its p as NA. popmix reported p = 1 and called
the fit good, for a model with infeasible weights that nothing had tested.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

DIR = Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_6pop"
REF = json.loads((DIR / "qpadm_saturated.json").read_text())


def _num(x: str) -> float:
    return math.nan if x == "NA" else float(x)


@pytest.fixture(scope="module")
def f2():
    return pm.read_admixtools2_f2(str(DIR / "f2_dir"))


@pytest.fixture(scope="module")
def result(f2):
    return pm.qpadm(f2, REF["left"], REF["right"], REF["target"])


def test_weights_match_admixtools2(result):
    np.testing.assert_allclose(result.weights, [_num(w["weight"]) for w in REF["weights"]], rtol=1e-6)
    np.testing.assert_allclose(result.weight_se, [_num(w["se"]) for w in REF["weights"]], rtol=1e-6)


def test_the_saturated_model_has_no_p_value(result):
    full = REF["rankdrop"][0]
    assert full["dof"] == 0 and full["p"] == "NA"
    assert result.dof == 0
    assert math.isnan(result.p)


def test_popdrop_matches_admixtools2(result):
    ours = {row.pattern: row for row in result.popdrop}
    for want in REF["popdrop"]:
        row = ours[want["pat"]]
        assert row.dof == want["dof"], want["pat"]
        if want["p"] == "NA":
            assert math.isnan(row.p), want["pat"]
        else:
            np.testing.assert_allclose(row.p, _num(want["p"]), rtol=1e-4, err_msg=want["pat"])


def _cli(f2, tmp_path: Path, fmt: str) -> str:
    d = tmp_path / "f2"
    f2.save(str(d))
    out = subprocess.run(["popmix", "qpadm", "--f2-dir", str(d), "--target", REF["target"],
                          "--left", ",".join(REF["left"]), "--right", ",".join(REF["right"]), "--format", fmt],
                         check=True, capture_output=True, text=True)
    return out.stdout


def test_cli_says_the_model_is_not_tested(f2, tmp_path):
    text = _cli(f2, tmp_path, "text")
    assert "Good model fit" not in text
    assert "=> Not tested: 0 degrees of freedom" in text
    assert "p=NA" in text


def test_cli_json_is_strict_json_with_a_null_p(f2, tmp_path):
    def no_nan(token: str):
        raise ValueError(f"not JSON: {token}")

    out = json.loads(_cli(f2, tmp_path, "json"), parse_constant=no_nan)
    assert out["p"] is None
    assert out["dof"] == 0
