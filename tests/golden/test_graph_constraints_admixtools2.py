"""Golden master tests: popmix graph constraints vs admixtools2.

The reference (tests/golden/admixtools2/v2.0.10/dg_6pop/graph_constraints.json,
from generate_graph_constraints.R) holds, for every 4th graph of each case
in all_graphs.json, admixtools2's summarize_numadmix, summarize_zerof4 and
the status (TRUE, FALSE or undecided) of each of several event-order
constraints of both types. All must match exactly.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import popmix as pm

REF_FILE = Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_6pop" / "graph_constraints.json"

pytestmark = pytest.mark.skipif(not REF_FILE.exists(), reason="graph_constraints.json not generated")


def _cases():
    if not REF_FILE.exists():
        return []
    cases = json.loads(REF_FILE.read_text())["cases"]
    return [pytest.param(c, id=f"{len(c['leaves'])}leaves_{c['nadmix']}admix") for c in cases]


def _edges(g: dict) -> list:
    return [tuple(e) for e in g["edges"]]


@pytest.mark.parametrize("case", _cases())
def test_summarize_numadmix(case):
    for g in case["graphs"]:
        assert pm.summarize_numadmix(_edges(g)) == g["numadmix"], g["index"]


@pytest.mark.parametrize("case", _cases())
def test_summarize_zerof4(case):
    for g in case["graphs"]:
        assert pm.summarize_zerof4(_edges(g)) == {tuple(q) for q in g["zerof4"]}, g["index"]


@pytest.mark.parametrize("case", _cases())
def test_event_order(case):
    for g in case["graphs"]:
        for e in g["events"]:
            earlier1, earlier2, later1, later2 = e["event"]
            for t in (1, 2):
                c = [{"earlier1": earlier1, "earlier2": earlier2, "later1": later1, "later2": later2, "type": t}]
                want = e[f"type{t}"]
                where = (g["index"], e["event"], t)
                assert pm.satisfies_eventorder(_edges(g), c) is (want is True), where
                assert pm.satisfies_eventorder(_edges(g), c, strict=False) is (want is not False), where


@pytest.mark.parametrize("case", _cases())
def test_constraints_from_the_summaries(case):
    """satisfies_numadmix and satisfies_(non)zerof4 follow from the summaries."""
    for g in case["graphs"]:
        edges, n = _edges(g), g["numadmix"]
        pop = case["leaves"][-1]
        assert pm.satisfies_numadmix(edges, [{"pop": pop, "min": n[pop], "max": None}])
        assert not pm.satisfies_numadmix(edges, [{"pop": pop, "min": n[pop] + 1, "max": None}])
        assert pm.satisfies_numadmix(edges, [{"pop": "absent", "min": 5, "max": 5}])
        zero = [tuple(q) for q in g["zerof4"]]
        assert pm.satisfies_zerof4(edges, zero)
        if zero:
            # the same statistics written another way round are still zero, and not nonzero
            a, b, c, d = zero[0]
            assert pm.satisfies_zerof4(edges, [(d, c, b, a)])
            assert not pm.satisfies_nonzerof4(edges, [(b, a, c, d)])
