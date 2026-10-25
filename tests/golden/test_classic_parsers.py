"""Golden master tests: popmix's parsers of classic ADMIXTOOLS output vs
admixtools2's.

tests/golden/admixtools_classic/v7.0.2 holds the output of ADMIXTOOLS 7.0.2
qpDstat, qp3Pop, qpAdm, qpF4ratio, qpGraph and qpfstats on dg_6pop
populations of AADR v62 (generate_classic.py) and, in parsed.json, what
admixtools2's parse_* functions make of each. popmix's must give the same
tables: values read from the files exactly, derived columns (se from a
statistic and its Z, p-values, z) within 1e-12.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import popmix as pm
import pytest

DIR = Path(__file__).parent / "admixtools_classic" / "v7.0.2"
pytestmark = pytest.mark.skipif(not (DIR / "parsed.json").exists(), reason="classic outputs not generated")


def _ref() -> dict:
    return json.loads((DIR / "parsed.json").read_text())


def _value(x):
    if isinstance(x, str):
        try:
            return float(x)
        except ValueError:
            return x
    return x


def _same_table(ours, theirs: list, label: str):
    assert len(ours) == len(theirs), label
    assert list(ours.columns) == list(theirs[0]), label
    for k, (row, want) in enumerate(zip(ours.to_dict("records"), theirs)):  # noqa: B905
        for col, w in want.items():
            o = row[col]
            w = w if isinstance(o, str) else _value(w)
            where = f"{label} row {k} {col}"
            if isinstance(w, float):
                if math.isnan(w):
                    assert math.isnan(o), where
                else:
                    assert o == pytest.approx(w, rel=1e-12, abs=1e-300), where
            else:
                assert o == w, where


def test_qpdstat():
    _same_table(pm.parse_qpdstat_output(DIR / "qpdstat.out"), _ref()["qpdstat"], "qpdstat")


def test_qp3pop():
    _same_table(pm.parse_qp3pop_output(DIR / "qp3pop.out"), _ref()["qp3pop"], "qp3pop")


def test_qpadm():
    ours = pm.parse_qpadm_output(DIR / "qpadm.out")
    _same_table(ours["weights"], _ref()["qpadm"]["weights"], "qpadm weights")
    _same_table(ours["popdrop"], _ref()["qpadm"]["popdrop"], "qpadm popdrop")


def test_qpf4ratio():
    _same_table(pm.parse_qpf4ratio_output(DIR / "qpf4ratio.out"), _ref()["qpf4ratio"], "qpf4ratio")


def test_qpgraph():
    ours = pm.parse_qpgraph_output(DIR / "qpgraph.ggg", DIR / "qpgraph.out")
    ref = _ref()["qpgraph"]
    assert ours["score"] == float(ref["score"])
    for k in ("edges", "f2", "f3"):
        _same_table(ours[k], ref[k], f"qpgraph {k}")
    assert (ours["f4"] is None) == (ref["f4"] is None)
    _same_table(pm.parse_qpgraph_output_edges(DIR / "qpgraph.ggg"), _ref()["qpgraph_edges_ggg"], "ggg edges")


def test_qpgraph_graphfile():
    labels = pm.parse_qpgraph_graphfile(DIR / "graph.txt", compat="admixtools2-2.0.10")
    _same_table(labels, _ref()["graphfile"], "graph file")
    pops = pm.parse_qpgraph_graphfile(DIR / "graph.txt")
    assert "Chimp.REF" in set(pops["to"]) and "Chimp" not in set(pops["to"])
    assert list(pm.parse_qpgraph_graphfile(DIR / "graph.txt", uselabels=True)["to"]) == list(labels["to"])


def test_qpfstats():
    ours = pm.parse_qpfstats_output(DIR / "qpfstats.out")
    ref = _ref()["qpfstats"]
    _same_table(ours["basis"], ref["basis"], "qpfstats basis")
    assert (ours["fstat"] is None) == (ref["fstat"] is None)
    assert (ours["jest"] is None) == (ref["jest"] is None)
