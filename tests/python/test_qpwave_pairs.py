"""The Rust qpwave_pairs result (popmix._core): its clades and best clade
come from its pairs (TASKS.md I4). The public ``pm.qpwave_pairs`` is the
admixtools2-shaped table (test_qpadm_models.py); which to keep is J2."""

import math
import os
from pathlib import Path

import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"
pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def _halves():
    """The French and the Yoruba each split in two populations."""
    inds, pops = [], []
    for line in Path(f"{SMALL}.ind").read_text().splitlines():
        sample, _, pop = line.split()
        if pop in ("French", "Yoruba"):
            inds.append(sample)
            pops.append(f"{pop}{len([p for p in pops if p.startswith(pop)]) % 2 + 1}")
        elif pop == "Simulation.1":
            inds.append(sample)
            pops.append(pop)
    return inds, pops


def test_clades_and_the_best_clade_follow_the_p_values():
    inds, pops = _halves()
    names = ["French1", "French2", "Simulation.1", "Yoruba1", "Yoruba2"]
    f2 = pm.extract_f2(SMALL, names, inds=inds, pops=pops)
    res = pm._core.qpwave_pairs(f2, ["French1", "French2", "Simulation.1"], ["Yoruba1", "Yoruba2"])


    assert res.n_pairs == 3
    pairs = res.pairs
    assert len(pairs) + res.n_skipped == 3
    for threshold in (0.0, 0.05, 0.5):
        want = [(p.pop_a, p.pop_b) for p in pairs if p.p_value > threshold]
        assert [(p.pop_a, p.pop_b) for p in res.clades(threshold)] == want
    candidates = [p for p in pairs if not math.isnan(p.p_value)]
    best = res.best_clade()
    if candidates:
        top = max(candidates, key=lambda p: p.p_value)
        assert (best.pop_a, best.pop_b, best.p_value) == (top.pop_a, top.pop_b, top.p_value)
    else:
        assert best is None
