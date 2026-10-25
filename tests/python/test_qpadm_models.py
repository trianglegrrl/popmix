"""Tests for the qpAdm model-set functions (popmix.qpadm_models).

Agreement with admixtools2 is in tests/golden/test_qpadm_models_admixtools2.py;
these check the interfaces and that the functions agree with qpadm and qpwave.
"""

import os

import pandas as pd
import popmix as pm
import pytest

PREFIX = "tests/golden/dates/family_chr20_22"
needs_data = pytest.mark.skipif(not os.path.exists(f"{PREFIX}.geno"), reason="test data not available")

TARGET = "Simulation.1"
RIGHT = ["Simulation", "Yoruba"]


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(PREFIX, populations=["Simulation", "French", "Yoruba", "Simulation.1"], auto_only=False)


def test_qpadm_models_splits_in_admixtools2_order():
    m = pm.qpadm_models(["A", "B", "C"])
    assert m.to_dict("records") == [{"left": ["A"], "right": ["B", "C"]}, {"left": ["B"], "right": ["A", "C"]},
                                    {"left": ["C"], "right": ["A", "B"]}]
    every = pm.qpadm_models(["A", "B", "C"], allpops=False, more_right=False)
    assert {"left": ["A", "B"], "right": ["C"]} in every.to_dict("records")
    assert {"left": ["A"], "right": ["B"]} in every.to_dict("records")
    assert all(len(r) > 0 for r in every.right)


def test_rotate_models_puts_the_rest_then_rightfix_on_the_right():
    m = pm.rotate_models(["A", "B", "C"], "T", "R")
    assert m.to_dict("records") == [
        {"left": ["A"], "right": ["B", "C", "R"], "target": "T"},
        {"left": ["B"], "right": ["A", "C", "R"], "target": "T"},
        {"left": ["C"], "right": ["A", "B", "R"], "target": "T"},
    ]
    # more fixed right populations allow more left ones
    assert pm.rotate_models(["A", "B", "C"], "T", ["R", "S"]).left.map(len).max() == 2


@needs_data
def test_qpadm_p_matches_qpadm(f2):
    full = pm.qpadm(f2, TARGET, ["French"], RIGHT)
    p = pm.qpadm_p(f2, ["French"], RIGHT, TARGET, weights=True)
    assert p["f4rank"] == 0
    assert p["dof"] == full.df
    assert p["chisq"] == pytest.approx(full.chi_sq, rel=1e-12)
    assert p["weights"] == {"French": 1.0}
    assert p["feasible"]
    assert "weights" not in pm.qpadm_p(f2, ["French"], RIGHT, TARGET)


@needs_data
def test_test_cladality_is_qpwave_rank_zero(f2):
    left = [TARGET, "French"]
    rank0 = next(t for t in pm.qpwave(f2, left, RIGHT).rank_tests if t.rank == 0)
    c = pm.test_cladality(f2, left, RIGHT)
    assert c["dof"] == rank0.df
    assert c["chisq"] == pytest.approx(rank0.chi_sq, rel=1e-12)


@needs_data
def test_qpwave_pairs_lists_each_pair_both_ways(f2):
    t = pm.qpwave_pairs(f2, ["French", TARGET], RIGHT)
    assert list(t.columns) == ["pop1", "pop2", "chisq", "p"]
    assert list(zip(t.pop1, t.pop2)) == [("French", TARGET), (TARGET, "French")]  # noqa: B905
    assert t.chisq[0] == t.chisq[1]
    assert t.chisq[0] == pytest.approx(pm.test_cladality(f2, ["French", TARGET], RIGHT)["chisq"], rel=1e-12)


@needs_data
def test_qpadm_multi_takes_any_table_of_models(f2):
    rows = [{"left": ["French"], "right": RIGHT, "target": TARGET}, {"left": [TARGET, "French"], "right": RIGHT}]
    as_dict = {"left": [r["left"] for r in rows], "right": [r["right"] for r in rows],
               "target": [r.get("target") for r in rows]}
    tables = [pm.qpadm_multi(f2, models, full_results=False) for models in (rows, as_dict, pd.DataFrame(as_dict))]
    for t in tables[1:]:
        pd.testing.assert_frame_equal(t, tables[0])
    # without a target the first left population is the target: the same model
    assert tables[0].chisq[0] == tables[0].chisq[1]
    full = pm.qpadm_multi(f2, rows)
    assert [r.sources for r in full] == [["French"], ["French"]]
    assert full[0].chi_sq == pytest.approx(tables[0].chisq[0], rel=1e-12)


def test_qpadm_multi_rejects_malformed_models():
    with pytest.raises(ValueError, match="equal length"):
        pm.qpadm_multi(None, {"left": [["A"], ["B"]], "right": [["C", "D"]]})
    with pytest.raises(ValueError, match="'left' and 'right'"):
        pm.qpadm_multi(None, [{"left": ["A"]}])


@needs_data
def test_qpadm_rotate_full_results_are_qpadm(f2):
    rows = pm.qpadm_rotate(f2, ["French", "Yoruba"], TARGET, "Simulation")
    full = pm.qpadm_rotate(f2, ["French", "Yoruba"], TARGET, "Simulation", full_results=True)
    assert rows.left.tolist() == [r.sources for r in full] == [["French"], ["Yoruba"]]
    assert rows.right.tolist() == [["Simulation", "Yoruba"], ["Simulation", "French"]]
    for (_, row), r in zip(rows.iterrows(), full):  # noqa: B905
        assert row.chisq == pytest.approx(r.chi_sq, rel=1e-12)
