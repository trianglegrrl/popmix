"""qpAdm's verdict weighs p, feasibility, the SEs and the simpler models (G1).

A model is "not rejected", never "a good fit", and only when p exceeds
alpha and every weight is in [0, 1]; the notes say when a weight is within
2 SE of 0 or a model with fewer sources is also not rejected. Rotations and
model sets say how many models they tested (Harney et al. 2021).
"""

from __future__ import annotations

import json
import subprocess
import warnings
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
POPS = ["French_1", "French_2", "French_3", "Simulation", "Simulation.1", "Yoruba_1", "Yoruba_2", "Yoruba_3"]
RIGHT = ["French_2", "Yoruba_2", "French_3", "Yoruba_3"]


@pytest.fixture(scope="module")
def family6(tmp_path_factory) -> str:
    from test_skill import example_data

    d = tmp_path_factory.mktemp("family6")
    example_data(d)
    return str(d / "family6")


@pytest.fixture(scope="module")
def f2(family6):
    return pm.extract_f2(family6, POPS)


def test_a_model_that_fits_is_not_rejected_and_has_no_notes(f2):
    a = pm.qpadm(f2, ["French_1", "Yoruba_1"], RIGHT, "Simulation.1").assess()
    assert a == pm.QpAdmAssessment(status="not rejected", alpha=0.05, notes=())


def test_a_small_p_rejects(f2):
    a = pm.qpadm(f2, ["French_1"], RIGHT, "Simulation.1").assess()
    assert a.status == "rejected"
    assert a.notes == ()


def test_alpha_is_the_level(f2):
    r = pm.qpadm(f2, ["French_2", "Yoruba_1"], ["Yoruba_2", "French_3", "Yoruba_3"], "French_1")
    assert r.p == pytest.approx(0.0243735, rel=1e-5)
    assert r.assess().status == "rejected"
    assert r.assess(alpha=0.01).status == "not rejected"


def test_a_weight_outside_0_1_is_infeasible_and_a_simpler_model_is_named(f2):
    a = pm.qpadm(f2, ["French_2", "French_3"], ["Yoruba_2", "Yoruba_1", "Yoruba_3"], "French_1").assess()
    assert a.status == "infeasible"
    assert a.notes == (
        "French_2: weight -0.223 is outside [0, 1]",
        "French_3: weight 1.223 is outside [0, 1]",
        "a model with fewer sources is also not rejected: French_3 alone (p = 0.857)",
    )


def test_a_weight_within_2_se_of_0_and_a_simpler_model_are_noted(f2):
    a = pm.qpadm(f2, ["Yoruba_1", "Simulation"], RIGHT, "Simulation.1").assess()
    assert a.status == "not rejected"
    assert a.notes == (
        "Yoruba_1: weight 0.073 is within 2 SE (0.238) of 0",
        "a model with fewer sources is also not rejected: Simulation alone (p = 0.405)",
    )


def test_the_simplest_of_several_simpler_models_is_named(f2):
    a = pm.qpadm(f2, ["French_1", "Yoruba_1", "Simulation"], RIGHT, "Simulation.1", constrained=True).assess()
    assert a.notes[-1] == ("4 models with fewer sources are also not rejected; the simplest: Simulation alone "
                           "(p = 0.346)")


def test_no_degrees_of_freedom_is_not_tested(f2):
    a = pm.qpadm(f2, ["French_1", "Yoruba_1"], ["French_2", "Yoruba_2"], "Simulation.1").assess()
    assert a.status == "not tested"


def test_cli_prints_the_assessment_not_a_good_fit(family6):
    out = subprocess.run(["popmix", "qpadm", "-i", family6, "--target", "Simulation.1", "--left",
                          "Yoruba_1,Simulation", "--right", ",".join(RIGHT)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert "Good model fit" not in out.stdout
    assert ("=> Not rejected at alpha 0.05, with every weight in [0, 1]. That does not make it the "
            "model: others may fit too.") in " ".join(out.stdout.split())
    assert "   - Yoruba_1: weight 0.073 is within 2 SE (0.238) of 0" in out.stdout


def test_cli_takes_alpha_and_writes_the_assessment_to_json(family6):
    out = subprocess.run(["popmix", "qpadm", "-i", family6, "--target", "French_1", "--left", "French_2,Yoruba_1",
                          "--right", "Yoruba_2,French_3,Yoruba_3", "--alpha", "0.01", "--format", "json"],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert json.loads(out.stdout)["assessment"] == {"status": "not rejected", "alpha": 0.01, "notes": []}


@pytest.mark.parametrize("call", ["rotate", "multi"])
def test_model_sets_report_how_many_models_they_tested(f2, call):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if call == "rotate":
            table = pm.qpadm_rotate(f2, ["French_1", "Yoruba_1", "Simulation"], "Simulation.1", RIGHT)
        else:
            table = pm.qpadm_multi(f2, [{"left": ["French_1"], "right": RIGHT, "target": "Simulation.1"},
                                        {"left": ["Yoruba_1"], "right": RIGHT, "target": "Simulation.1"}],
                                   full_results=False)
    n = len(table)
    assert table.attrs["n_models"] == n
    messages = [str(w.message) for w in caught if issubclass(w.category, pm.PopmixWarning)]
    assert messages == [f"{n} qpAdm models tested: their p-values are not corrected for multiple testing, so "
                        "report how many models were tried with any model chosen (Harney et al. 2021)"]


def test_qpadm_multi_returns_result_objects_that_assess(f2):
    results = pm.qpadm_multi(f2, [{"left": ["French_1", "Yoruba_1"], "right": RIGHT, "target": "Simulation.1"}])
    assert isinstance(results[0], pm.QpAdmResult)
    assert results[0].assess().status == "not rejected"
