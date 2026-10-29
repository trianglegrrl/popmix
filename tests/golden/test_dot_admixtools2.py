"""Golden master tests: DOT output vs admixtools2 write_dot / parse_dot.

dot.json (generate_dot.R) holds admixtools2's write_dot() text for the
dg_6pop qpgraph fits (fitted edges, with a title) and for the bare graphs
(another page size), and parse_dot() of each. popmix with
compat="admixtools2-2.0.10" must write the same text; its own output keeps
the names, and parses back to the graph's edges.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import popmix as pm

from .reference_paths import missing

REF = Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_6pop" / "dot.json"
pytestmark = pytest.mark.skipif(missing(REF), reason="dot.json not generated")


def _graphs():
    if missing(REF):
        return []
    return [pytest.param(g, id=g["name"]) for g in json.loads(REF.read_text(encoding="utf-8"))["graphs"]]


@pytest.mark.parametrize("g", _graphs())
def test_compat_text_is_admixtools2s(g):
    fit = [(e["from"], e["to"], e["type"], e["weight"]) for e in g["fit_edges"]]
    assert pm.write_dot(fit, title=g["name"], compat="admixtools2-2.0.10") == g["fitted"]["text"] + "\n"
    bare = [(e["from"], e["to"]) for e in g["fit_edges"]]
    assert pm.write_dot(bare, size=(5, 4.25), compat="admixtools2-2.0.10") == g["bare"]["text"] + "\n"


@pytest.mark.parametrize("g", _graphs())
def test_parse_dot(g):
    assert [list(e) for e in pm.parse_dot(g["fitted"]["text"])] == g["fitted"]["parsed"]
    edges = [(e["from"], e["to"]) for e in g["fit_edges"]]
    assert pm.parse_dot(pm.write_dot(edges)) == edges


def test_qpgraph_result_and_file(tmp_path):
    g = _graphs()[-1].values[0]
    fit = [(e["from"], e["to"], e["type"], e["weight"]) for e in g["fit_edges"]]
    text = pm.write_dot(fit, path=tmp_path / "g.dot")
    assert (tmp_path / "g.dot").read_text(encoding="utf-8") == text
    assert '"Han.DG" [shape' not in text and '"n1" [shape = point];' in text
    assert pm.parse_dot(tmp_path / "g.dot") == [(a, b) for a, b, *_ in fit]
    with pytest.raises(ValueError, match="compat"):
        pm.write_dot(fit, compat="admixtools2-9")
