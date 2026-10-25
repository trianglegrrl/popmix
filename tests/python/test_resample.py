"""Tests for resampling (popmix.resample).

Agreement with admixtools2 is in tests/golden/test_resample_admixtools2.py;
these check the interfaces on the small test data.
"""

import os
from pathlib import Path

import pytest

import popmix as pm

PREFIX = "tests/golden/dates/family_chr20_22"
pytestmark = pytest.mark.skipif(not os.path.exists(f"{PREFIX}.geno"), reason="test data not available")

RIGHT = ["Simulation", "Yoruba"]


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(PREFIX, populations=["Simulation", "French", "Yoruba", "Simulation.1"], auto_only=False)


def test_loo_list_leaves_out_each_block(f2):
    loo = pm.loo_list(f2)
    assert len(loo) == f2.n_blocks
    for i, b in enumerate(loo):
        assert b.block_lengths == f2.block_lengths[:i] + f2.block_lengths[i + 1:]


def test_boo_list_is_reproducible_from_its_seed(f2):
    a, b, c = (pm.boo_list(f2, nboot=3, seed=s) for s in (5, 5, 6))
    assert a["sel"] == b["sel"] != c["sel"]
    assert all(len(s) == f2.n_blocks for s in a["sel"])
    assert all(t.n_blocks == f2.n_blocks - len(set(s)) for t, s in zip(a["test"], a["sel"]))  # noqa: B905


def test_qpadm_resample_snps_jackknife(f2):
    full = pm.qpadm(f2, "Simulation.1", ["French"], RIGHT)
    t = pm.qpadm_resample_snps(f2, "Simulation.1", ["French"], RIGHT)
    assert list(t.id) == list(range(1, f2.n_blocks + 1))
    assert t.error.isna().all()
    assert all(w == [1.0] for w in t.weights)
    # leaving out one block barely moves the statistic
    assert t.chisq.median() == pytest.approx(full.chi_sq, rel=0.2)


def test_bootstrap_count_and_seed(f2):
    a = pm.qpwave_resample_snps(f2, ["Simulation.1", "French"], RIGHT, boot=4, boot_seed=1)
    b = pm.qpwave_resample_snps(f2, ["Simulation.1", "French"], RIGHT, boot=4, boot_seed=1)
    assert len(a) == 4
    assert list(a.sel) == list(b.sel)
    assert [r[0]["chisq"] for r in a.rankdrop] == [r[0]["chisq"] for r in b.rankdrop]
    assert [r[0]["f4rank"] for r in a.rankdrop] == [0] * 4


def test_failures_are_recorded_per_resample(f2):
    t = pm.qpadm_resample_snps(f2, "Simulation.1", ["French"], ["Simulation", "Nowhere"], boot=2, boot_seed=1)
    assert len(t) == 2
    assert all(isinstance(e, Exception) for e in t.error)


def test_fstat_resample_snps_has_one_row_per_statistic_and_resample(f2):
    t = pm.qpdstat_resample_snps(f2, "Simulation", "French", ["Yoruba", "Simulation.1"], "Simulation.1", boot=3,
                                 boot_seed=2)
    assert len(t) == 3 * 2
    assert list(t.columns[:2]) == ["id", "sel"]
    assert list(t.id) == [1, 1, 2, 2, 3, 3]


def test_qpgraph_resample_multi_scores_held_out_blocks(f2):
    tree = [("R", "Simulation"), ("R", "n1"), ("n1", "French"), ("n1", "n2"), ("n2", "Yoruba"),
            ("n2", "Simulation.1")]
    fits = pm.qpgraph_resample_multi(f2, [tree, tree], nboot=3, boot_seed=3)
    assert len(fits) == 2 and all(len(f) == 3 for f in fits)
    assert list(fits[0].score_test) == list(fits[1].score_test)
    assert all(s > 0 for s in fits[0].score_test)


def test_resample_inds_leaves_out_each_individual_of_a_larger_population():
    ind = [line.split() for line in Path(f"{PREFIX}.ind").read_text().splitlines()]
    def first(pop, n):
        return [r for r in ind if r[2] == pop][:n]

    rows = first("French", 3) + first("Simulation.1", 2) + first("Simulation", 1) + first("Yoruba", 2)
    inds, pops = [r[0] for r in rows], [r[2] for r in rows]
    t = pm.qpadm_resample_inds(PREFIX, inds, pops, "Simulation.1", ["French"], RIGHT, f2_args={"auto_only": False})
    singles = [i for i, p in zip(inds, pops) if pops.count(p) == 1]  # noqa: B905
    assert list(t.ind) == [i for i in inds if i not in singles]
    assert t.error.isna().all()
    with pytest.raises(ValueError, match="same length"):
        pm.qpadm_resample_inds(PREFIX, inds, pops[:-1], "Simulation.1", ["French"], RIGHT)
