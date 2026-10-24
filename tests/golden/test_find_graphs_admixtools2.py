"""Graph search against admixtools2 2.0.10 generate_all_graphs() and qpgraph().

all_graphs.json (generate_all_graphs.R) holds, for several (leaves, nadmix)
cases on the dg_6pop f2 blocks, every graph generate_all_graphs() returns in
its order, each with its graph_hash() key and qpgraph() score. popmix must
generate the same graphs in the same order and, at admixtools2's admixture
weights, fit each to admixtools2's score; its searches must find a graph at
least as good as admixtools2's best.
"""

import json
import os
from pathlib import Path

import pytest

import popmix as pm

REF_DIR = Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_6pop"
REF_FILE = REF_DIR / "all_graphs.json"

pytestmark = pytest.mark.skipif(not REF_FILE.exists(), reason="all_graphs.json not generated")


def _reference() -> dict:
    return json.loads(REF_FILE.read_text())


def _cases():
    if not REF_FILE.exists():
        return []
    return [pytest.param(i, id=f"{len(c['leaves'])}leaves_{c['nadmix']}admix")
            for i, c in enumerate(_reference()["cases"])]


def _edges(g: dict) -> list:
    return [tuple(e) for e in g["edges"]]


@pytest.fixture(scope="module")
def f2_blocks():
    f2_ref = json.loads((REF_DIR / "f2.json").read_text())
    if not os.path.exists(f"{f2_ref['prefix']}.geno"):
        pytest.skip(f"genotype data not available: {f2_ref['prefix']}")
    assert f2_ref["pops"] == _reference()["pops"]
    p = f2_ref["params"]
    return pm.extract_f2(f2_ref["prefix"], f2_ref["pops"], auto_only=p["auto_only"], block_size=p["blgsize"],
                         maxmiss=p["maxmiss"], poly_only=p["poly_only"] == "f2")


# At fixed admixture weights the fit has no optimizer in it, so it is exact.
FIXED_RTOL = 1e-10


def _close(ours: float, theirs: float, rtol: float) -> bool:
    return abs(ours - theirs) <= rtol * abs(theirs) + 1e-15


@pytest.mark.parametrize("i", _cases())
def test_same_graphs_in_the_same_order(i):
    case = _reference()["cases"][i]
    ours = pm.all_graphs(case["leaves"], case["nadmix"])
    assert len(ours) == len(case["graphs"])
    ours_keys = [pm.graph_key(g) for g in ours]
    theirs = [g["key"] for g in case["graphs"]]
    diff = [k for k, (a, b) in enumerate(zip(ours_keys, theirs)) if a != b]
    assert not diff, f"{len(diff)} graphs differ, first at {diff[0]}"


@pytest.mark.parametrize("i", _cases())
def test_graph_key_of_admixtools2_graphs(i):
    case = _reference()["cases"][i]
    for g in case["graphs"]:
        assert pm.graph_key(_edges(g)) == g["key"]


@pytest.mark.slow
@pytest.mark.parametrize("i", _cases())
def test_qpgraph_at_admixtools2_weights(f2_blocks, i):
    """Every graph, fitted at admixtools2's best weights, scores what admixtools2 scores."""
    case = _reference()["cases"][i]
    bad = []
    for k, g in enumerate(case["graphs"]):
        weights = g["weights"] or None
        ours = pm.qpgraph(f2_blocks, _edges(g), admix_weights=weights).score
        if not _close(ours, g["score"], FIXED_RTOL):
            bad.append(f"graph {k}: {ours!r} != {g['score']!r}")
    assert not bad, f"{len(bad)} of {len(case['graphs'])} differ:\n" + "\n".join(bad[:10])


@pytest.mark.slow
@pytest.mark.parametrize("i", [c for c in _cases() if c.id.endswith("_0admix")])
def test_qpgraph_tree_scores(f2_blocks, i):
    """Trees have no weights to optimize: qpgraph() itself agrees exactly."""
    case = _reference()["cases"][i]
    for g in case["graphs"]:
        assert _close(pm.qpgraph(f2_blocks, _edges(g)).score, g["score"], FIXED_RTOL)


@pytest.mark.slow
@pytest.mark.xfail(strict=True, reason="weight optimizers find different local optima on some graphs; "
                                       "DEVIATIONS.md, qpGraph weight optimizer")
@pytest.mark.parametrize("i", [c for c in _cases() if not c.id.endswith("_0admix")])
def test_qpgraph_optimized_score_of_every_graph(f2_blocks, i):
    case = _reference()["cases"][i]
    worse = [k for k, g in enumerate(case["graphs"])
             if pm.qpgraph(f2_blocks, _edges(g)).score > g["score"] * (1 + 1e-6) + 1e-12]
    assert not worse, f"popmix's optimum is worse than admixtools2's for {len(worse)} graphs"


@pytest.mark.slow
@pytest.mark.parametrize("i", _cases())
def test_exhaustive_search(f2_blocks, i):
    """find_graphs(method="exhaustive") fits every graph, and its best is at least admixtools2's best."""
    case = _reference()["cases"][i]
    graphs, n_evaluated = pm.find_graphs(f2_blocks, case["leaves"], nadmix=case["nadmix"], method="exhaustive")
    assert n_evaluated == len(case["graphs"])
    ref = {tuple(g["key"]): g["score"] for g in case["graphs"]}
    ours = {tuple(pm.graph_key(g.edges)): g.score for g in graphs}
    assert ours.keys() == ref.keys()
    scores = [g.score for g in graphs]
    assert scores == sorted(scores)
    best = min(ref.values())
    if case["nadmix"] == 0:
        assert all(_close(ours[k], ref[k], FIXED_RTOL) for k in ref)
        assert _close(graphs[0].score, best, FIXED_RTOL)
    else:
        assert graphs[0].score <= best * (1 + 1e-6)


@pytest.mark.slow
def test_stepwise_search_finds_the_exhaustive_optimum(f2_blocks):
    """5 populations, 1 admixture: stepwise reaches the best score of all 1485 graphs."""
    case = next(c for c in _reference()["cases"] if len(c["leaves"]) == 5 and c["nadmix"] == 1)
    exhaustive, _ = pm.find_graphs(f2_blocks, case["leaves"], nadmix=1, method="exhaustive")
    graphs, n_evaluated = pm.find_graphs(f2_blocks, case["leaves"], nadmix=1)
    assert n_evaluated < len(case["graphs"])
    assert _close(graphs[0].score, exhaustive[0].score, FIXED_RTOL)
    assert graphs[0].score <= min(g["score"] for g in case["graphs"]) * (1 + 1e-6)


@pytest.mark.slow
def test_outgroup_restricts_to_graphs_with_the_outgroup_first(f2_blocks):
    case = next(c for c in _reference()["cases"] if len(c["leaves"]) == 5 and c["nadmix"] == 0)
    outgroup = case["leaves"][0]
    graphs, _ = pm.find_graphs(f2_blocks, case["leaves"], outgroup=outgroup, method="exhaustive")
    for g in graphs:
        root_children = [b for a, b in g.edges if a == g.edges[0][0]]
        assert outgroup in root_children
    # a tree on 5 leaves with a fixed first split: the 15 rooted trees of the other 4
    assert len(graphs) == 15
    with_outgroup = [g["score"] for g in case["graphs"]
                     if outgroup in [b for a, b in g["edges"] if a == g["edges"][0][0]]]
    assert len(with_outgroup) == 15
    assert _close(graphs[0].score, min(with_outgroup), FIXED_RTOL)
