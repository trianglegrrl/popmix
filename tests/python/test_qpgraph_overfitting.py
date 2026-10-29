"""qpGraph says when a fit overfits or its weights are not identifiable (G9).

A graph's score falls with every parameter added, so a good score is not
evidence: popmix reports the parameters against the f3 statistics fitted,
an AIC-style penalized score, an out-of-sample score (fit on alternate
jackknife blocks, scored on the others) by default on the CLI, and warns
when starts equally good reach different admixture weights.
"""

from __future__ import annotations

import json
import subprocess
import warnings
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
DG6 = str(ROOT / "tests" / "golden" / "admixtools2" / "v2.0.10" / "dg_6pop" / "f2_dir")
REF = json.loads((ROOT / "tests/golden/admixtools2/v2.0.10/dg_6pop/qpgraph.json").read_text(encoding="utf-8"))
ADMIX1 = [tuple(e) for e in next(g for g in REF["graphs"] if g["name"] == "admix1")["edges"]]
# two admixtures on Mbuti, French, Han, Papuan whose weights trade off: every
# start reaches the same score with a different weight
UNIDENTIFIED = [("R", "n1"), ("n1", "Mbuti.DG"), ("n2", "Han.DG"), ("n2", "n3"), ("n3", "Papuan.DG"), ("n1", "n4"),
                ("n3", "n4"), ("R", "n5"), ("n5", "n2"), ("n4", "n6"), ("n6", "French.DG"), ("n5", "n6")]


def caught(f) -> list[str]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        f()
    return [str(x.message) for x in w if issubclass(x.category, pm.PopmixWarning)]


@pytest.fixture(scope="module")
def dg6():
    return pm.read_admixtools2_f2(DG6)


def test_the_fit_counts_its_parameters(dg6):
    r = pm.qpgraph(dg6, ADMIX1)
    n_leaves = len(r.pops)
    assert r.n_f3 == n_leaves * (n_leaves - 1) // 2
    drift = sum(e.kind != "admix" for e in r.edges)
    admix = sum(e.kind == "admix" for e in r.edges) // 2
    assert r.n_params == drift + admix
    assert r.dof == r.n_f3 - r.n_params
    assert r.penalized_score == pytest.approx(r.score + 2 * r.n_params, rel=1e-15)


def test_the_out_of_sample_score_fits_one_half_of_the_blocks_and_scores_the_other(dg6):
    n = dg6.n_blocks
    even, odd = dg6.select_blocks(list(range(0, n, 2))), dg6.select_blocks(list(range(1, n, 2)))
    a = pm.qpgraph(even, ADMIX1, f2_blocks_test=odd).score_test
    b = pm.qpgraph(odd, ADMIX1, f2_blocks_test=even).score_test
    assert pm.qpgraph_out_of_sample(dg6, ADMIX1) == pytest.approx((a + b) / 2, rel=1e-12)


def test_unidentifiable_weights_warn(dg6):
    messages = caught(lambda: pm.qpgraph(dg6, UNIDENTIFIED))
    assert len(messages) == 1
    assert messages[0].startswith("qpgraph: 10 starts reach the best score (within 0.01%) with weights from ")
    assert messages[0].endswith(": the admixture weights are not identifiable from these f3 statistics")


def test_identifiable_weights_do_not_warn(dg6):
    assert caught(lambda: pm.qpgraph(dg6, ADMIX1)) == []


def test_the_cli_reports_parameters_and_the_out_of_sample_score(dg6, tmp_path):
    f2_dir = str(tmp_path / "f2")
    dg6.save(f2_dir)
    graph = tmp_path / "graph.txt"
    graph.write_text("".join(f"{a}\t{b}\n" for a, b in ADMIX1), encoding="utf-8")
    r = pm.qpgraph(dg6, ADMIX1)
    oos = pm.qpgraph_out_of_sample(dg6, ADMIX1)
    out = subprocess.run(["popmix", "qpgraph", "--f2-dir", f2_dir, "--graph", str(graph)], capture_output=True,
                         text=True)
    assert out.returncode == 0, out.stderr
    assert (f"Parameters: {r.n_params} for {r.n_f3} f3 statistics ({r.dof} dof); penalized score "
            f"{r.penalized_score:.4f}") in out.stdout
    assert f"Out-of-sample score (alternate blocks): {oos:.4f}" in out.stdout
    js = json.loads(subprocess.run(["popmix", "qpgraph", "--f2-dir", f2_dir, "--graph", str(graph), "--format",
                                    "json"], capture_output=True, text=True).stdout)
    assert js["n_params"] == r.n_params
    assert js["score_out_of_sample"] == pytest.approx(oos, rel=1e-12)
    quick = subprocess.run(["popmix", "qpgraph", "--f2-dir", f2_dir, "--graph", str(graph), "--no-out-of-sample"],
                           capture_output=True, text=True).stdout
    assert "Out-of-sample" not in quick


def test_find_graphs_reports_the_top_graphs_out_of_sample(dg6, tmp_path):
    f2_dir = str(tmp_path / "f2")
    dg6.save(f2_dir)
    out = subprocess.run(["popmix", "find-graphs", "--f2-dir", f2_dir, "--populations",
                          "Mbuti.DG,French.DG,Han.DG,Papuan.DG", "--nadmix", "1", "--top", "2", "--format", "json"],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    graphs = json.loads(out.stdout)["graphs"]
    for g in graphs:
        edges = [tuple(e) for e in g["edges"]]
        assert g["score_out_of_sample"] == pytest.approx(pm.qpgraph_out_of_sample(dg6, edges), rel=1e-12)
        assert g["n_params"] == pm.qpgraph(dg6, edges).n_params
