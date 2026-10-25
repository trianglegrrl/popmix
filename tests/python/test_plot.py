"""Tests for popmix.plot's graph, qpGraph fit and ROLLOFF plots (the others
are tested with their analyses)."""

import os

import popmix as pm
import pytest

plot = pytest.importorskip("popmix.plot")
pytest.importorskip("matplotlib")

PREFIX = "tests/golden/dates/family_chr20_22"
needs_data = pytest.mark.skipif(not os.path.exists(f"{PREFIX}.geno"), reason="test data not available")
GRAPH = [("R", "Yoruba"), ("R", "n1"), ("n1", "French"), ("n1", "n2"), ("n2", "Simulation"), ("n2", "a"),
         ("n1", "a"), ("a", "Simulation.1")]


def test_layout_puts_children_below_parents():
    pos = plot._graph_layout(GRAPH)
    assert all(pos[a][1] > pos[b][1] for a, b in GRAPH)
    leaves = ["Yoruba", "French", "Simulation", "Simulation.1"]
    assert len({pos[v][0] for v in leaves}) == len(leaves)


def test_graph_from_edges(tmp_path):
    ax = plot.graph(GRAPH, path=tmp_path / "g.png")
    assert (tmp_path / "g.png").stat().st_size > 0
    labels = {t.get_text() for t in ax.texts}
    assert {"Yoruba", "French", "Simulation", "Simulation.1"} <= labels


@needs_data
def test_graph_and_fit_of_a_qpgraph_result(tmp_path):
    f2 = pm.extract_f2(PREFIX, ["Yoruba", "French", "Simulation", "Simulation.1"], auto_only=False)
    fit = pm.qpgraph(f2, GRAPH)
    ax = plot.graph(fit, path=tmp_path / "g.png")
    texts = [t.get_text() for t in ax.texts]
    assert sum(t.endswith("%") for t in texts) == 2
    ax = plot.qpgraph_fit(fit, path=tmp_path / "f.png")
    assert len(ax.get_yticklabels()) == len(fit.f3)


@needs_data
def test_rolloff_curve(tmp_path):
    r = pm.rolloff(PREFIX, "Simulation.1", ["French", "Yoruba"])
    plot.rolloff_curve(r, path=tmp_path / "r.png")
    assert (tmp_path / "r.png").stat().st_size > 0
