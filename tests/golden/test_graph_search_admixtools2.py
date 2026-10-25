"""Golden master tests: the evolutionary graph search vs admixtools2 find_graphs().

find_graphs.json (generate_find_graphs.R) holds admixtools2 find_graphs()
on the dg_6pop f2 blocks, Chimp.REF as outgroup, default settings, for 5
populations with 1 admixture and 4 with 2, from 8 seeds: the best score of
each run and the number of graphs it fitted. The two searches use different
random numbers, so runs are compared as a group: over 30 seeds popmix must
reach the best score of admixtools2's runs (and the exhaustive optimum of
all_graphs.json), and with one admixture do so about as often.
"""

from __future__ import annotations

import json
import os
import statistics
from pathlib import Path

import popmix as pm
import pytest

REF_DIR = Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_6pop"
REF_FILE = REF_DIR / "find_graphs.json"
SEEDS = range(1, 31)

pytestmark = [pytest.mark.slow, pytest.mark.skipif(not REF_FILE.exists(), reason="find_graphs.json not generated")]


def _cases():
    if not REF_FILE.exists():
        return []
    return [pytest.param(c, id=f"{len(c['leaves'])}leaves_{c['nadmix']}admix")
            for c in json.loads(REF_FILE.read_text())["cases"]]


@pytest.fixture(scope="module")
def f2_blocks():
    f2_ref = json.loads((REF_DIR / "f2.json").read_text())
    if not os.path.exists(f"{f2_ref['prefix']}.geno"):
        pytest.skip(f"genotype data not available: {f2_ref['prefix']}")
    p = f2_ref["params"]
    return pm.extract_f2(f2_ref["prefix"], f2_ref["pops"], auto_only=p["auto_only"], block_size=p["blgsize"],
                         maxmiss=p["maxmiss"], poly_only=p["poly_only"] == "f2")


_runs: dict = {}


def _search(f2_blocks, case) -> list[tuple[float, int]]:
    key = (tuple(case["leaves"]), case["nadmix"])
    if key not in _runs:
        _runs[key] = []
        for s in SEEDS:
            graphs, n = pm.find_graphs(f2_blocks, case["leaves"], nadmix=case["nadmix"], outgroup=case["leaves"][0],
                                       method="evolutionary", seed=s)
            _runs[key].append((graphs[0].score, n))
    return _runs[key]


def _exhaustive_best(case) -> float:
    ref = json.loads((REF_DIR / "all_graphs.json").read_text())
    (c,) = [c for c in ref["cases"] if c["leaves"] == case["leaves"] and c["nadmix"] == case["nadmix"]]
    return min(g["score"] for g in c["graphs"])


@pytest.mark.parametrize("case", _cases())
def test_reaches_admixtools2s_best(f2_blocks, case):
    ours = min(s for s, _ in _search(f2_blocks, case))
    theirs = min(float(r["best"]) for r in case["runs"])
    assert ours <= theirs * (1 + 1e-6)
    assert ours <= _exhaustive_best(case) * (1 + 1e-2)


@pytest.mark.parametrize("case", [c for c in _cases() if c.values[0]["nadmix"] == 1])
def test_reaches_the_optimum_as_often(f2_blocks, case):
    best = _exhaustive_best(case)
    ours = [abs(s - best) <= 1e-6 * best for s, _ in _search(f2_blocks, case)]
    theirs = [abs(float(r["best"]) - best) <= 1e-6 * best for r in case["runs"]]
    assert sum(ours) / len(ours) >= 0.4, f"popmix {sum(ours)}/{len(ours)}, admixtools2 {sum(theirs)}/{len(theirs)}"


@pytest.mark.parametrize("case", _cases())
def test_fits_a_similar_number_of_graphs(f2_blocks, case):
    ours = statistics.median(n for _, n in _search(f2_blocks, case))
    theirs = statistics.median(r["n_fitted"] for r in case["runs"])
    assert theirs / 2 <= ours <= theirs * 2, (ours, theirs)
