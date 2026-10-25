"""Tests for the evolutionary graph search and its graph operations.

Agreement with admixtools2 is in tests/golden/test_graph_search_admixtools2.py
and tests/golden/test_graph_constraints_admixtools2.py.
"""

import os
import random

import popmix as pm
import popmix.graph_search as gs
import pytest

PREFIX = "tests/golden/dates/family_chr20_22"
POPS = ["Yoruba", "French", "Simulation", "Simulation.1"]


def _graphs(n: int, nadmix: int, seed: int = 1):
    rng = random.Random(seed)
    leaves = ["O", "A", "B", "C", "D"]
    return [gs.random_admixturegraph(leaves, nadmix, rng, outgroup="O") for _ in range(n)], leaves


@pytest.mark.parametrize("nadmix", [0, 1, 2])
def test_random_graphs_are_valid(nadmix):
    graphs, leaves = _graphs(30, nadmix)
    for g in graphs:
        assert gs._ok(g, leaves)
        assert gs._numadmix(g) == nadmix
        assert gs._outpop(g) == "O"


@pytest.mark.parametrize("name", list(gs.MUTATIONS))
def test_mutations_keep_leaves_and_outgroup(name):
    graphs, leaves = _graphs(25, 1, seed=2)
    rng = random.Random(3)
    made = 0
    for g in graphs:
        out = gs.MUTATIONS[name](g, gs._Names(g), rng, fix_outgroup=True)
        if out is None:
            continue
        made += 1
        assert gs._ok(out, leaves), (name, g, out)
        if name != "place_root_random":
            assert gs._outpop(out) == "O", (name, out)
    assert made >= 10, f"{name} made only {made} graphs"


def test_negative_edge_repairs_give_valid_graphs():
    graphs, leaves = _graphs(25, 2, seed=4)
    rng = random.Random(5)
    for g in graphs:
        for a in gs._admix_nodes(g):
            for p in gs._parents(g, a):
                for f in (gs.rearrange_negadmix3, gs.replace_admix_with_random):
                    out = f(g, gs._Names(g), rng, p, a)
                    assert out is None or gs._ok(out, leaves)
        for frm, to in g:
            out = gs.rearrange_negdrift(g, gs._Names(g), rng, frm, to)
            assert out is None or gs._ok(out, leaves)
        for out in gs.swap_negdrift(g, gs._Names(g), rng, g[:4]):
            assert gs._ok(out, leaves)


def test_delete_and_insert_admix_change_the_count_by_one():
    graphs, leaves = _graphs(20, 1, seed=6)
    rng = random.Random(7)
    for g in graphs:
        names = gs._Names(g)
        more = gs._insert_admix(g, names, rng)
        assert more is not None and gs._numadmix(more) == 2 and gs._ok(more, leaves)
        (a,) = gs._admix_nodes(g)
        less = gs._delete_admix(g, names, gs._parents(g, a)[0], a)
        assert less is not None and gs._numadmix(less) == 0 and gs._ok(less, leaves)


def test_base_first_puts_the_base_first():
    graphs, _ = _graphs(5, 1)
    for g in graphs:
        e = gs.base_first(g, "C")
        assert sorted(e) == sorted(g)
        assert gs._leaves(e)[0] == "C"


needs_data = pytest.mark.skipif(not os.path.exists(f"{PREFIX}.geno"), reason="test data not available")


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(PREFIX, populations=POPS, auto_only=False)


@needs_data
def test_evolutionary_search_is_reproducible(f2):
    a, na = pm.find_graphs(f2, POPS, nadmix=1, outgroup="Yoruba", method="evolutionary", seed=3, stop_gen=6)
    b, nb = pm.find_graphs(f2, POPS, nadmix=1, outgroup="Yoruba", method="evolutionary", seed=3, stop_gen=6)
    assert na == nb == len(a)
    assert [(g.edges, g.score) for g in a] == [(g.edges, g.score) for g in b]
    assert [g.score for g in a] == sorted(g.score for g in a)
    assert len({tuple(pm.graph_key(g.edges)) for g in a}) == len(a)
    for g in a:
        assert gs._outpop(g.edges) == "Yoruba"
        assert gs._leaves(g.edges)[0] == "Yoruba"


@needs_data
def test_evolutionary_search_reaches_the_exhaustive_optimum(f2):
    exhaustive, _ = pm.find_graphs(f2, POPS, nadmix=1, outgroup="Yoruba", method="exhaustive")
    best = min(min(g.score for g in pm.find_graphs(f2, POPS, nadmix=1, outgroup="Yoruba", method="evolutionary",
                                                   seed=s)[0]) for s in (1, 2, 3))
    assert best <= exhaustive[0].score * (1 + 1e-6)


@needs_data
def test_admix_constraints_hold_once_met(f2):
    cons = [{"pop": "Simulation.1", "min": 1, "max": None}]
    graphs, _ = pm.find_graphs(f2, POPS, nadmix=1, outgroup="Yoruba", method="evolutionary", seed=2,
                               admix_constraints=cons, stop_gen=10)
    met = sorted((g for g in graphs if pm.satisfies_numadmix(g.edges, cons)), key=lambda g: g.generation)
    assert met, "no graph met the constraint"
    first = met[0].generation
    assert all(pm.satisfies_numadmix(g.edges, cons) for g in graphs if g.generation > first)


@needs_data
def test_initgraph_and_max_admix(f2):
    tree = [("R", "Yoruba"), ("R", "a"), ("a", "French"), ("a", "b"), ("b", "Simulation"), ("b", "Simulation.1")]
    graphs, _ = pm.find_graphs(f2, POPS, method="evolutionary", initgraph=tree, max_admix=1, seed=1,
                               stop_gen=12, plusminus_generations=2)
    assert {gs._numadmix(g.edges) for g in graphs} <= {0, 1}
    with pytest.raises(ValueError, match="leaves"):
        pm.find_graphs(f2, POPS, method="evolutionary", initgraph=tree[:-1] + [("b", "Other")])





@pytest.mark.slow
def test_msprime_script_runs(tmp_path):
    """A script from msprime_sim runs and writes genotypes popmix reads with
    the graph's populations (the simulation itself is unseeded, as
    admixtools2 writes it)."""
    pytest.importorskip("msprime")
    pytest.importorskip("multiprocess")
    edges = [("R", "O"), ("R", "n1"), ("n1", "A"), ("n1", "n2"), ("n2", "B"), ("n2", "C")]
    path = pm.msprime_sim(edges, outpref=tmp_path / "sim", nsnps=300, ind_per_pop=4, neff=200, mutation_rate=1e-2,
                          numcores=2, fix_leaf=True, run=True)
    assert os.path.exists(path)
    assert len((tmp_path / "sim.geno").read_text().splitlines()) == 300
    assert len((tmp_path / "sim.ind").read_text().splitlines()) == 16
    # the script spreads the SNPs over 50 Morgans: 10 blocks
    f2 = pm.extract_f2(str(tmp_path / "sim"), ["O", "A", "B", "C"], auto_only=False, block_size=5)
    assert f2.n_blocks >= 5
    assert all(f2.f2(a, b)["estimate"] > 0 for a, b in [("O", "A"), ("A", "B"), ("B", "C")])
