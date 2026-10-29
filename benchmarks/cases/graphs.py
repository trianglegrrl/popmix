"""Benchmark cases: qpGraph and find_graphs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

from ._common import Case, _r, subset
from .samples import SPREAD_SEEDS, _seed_runs

GRAPH = [("R", "Chimp"), ("R", "n1"), ("n1", "Mbuti"), ("n1", "n2"), ("n2", "Han"), ("n2", "n3"),
         ("n3", "n4"), ("n3", "n5"), ("n4", "Turkey_N"), ("n4", "a"), ("n5", "Russia_Samara_EBA_Yamnaya"),
         ("n5", "a"), ("a", "French")]


GRAPH_POPS = ["Chimp", "Mbuti", "Han", "Turkey_N", "Russia_Samara_EBA_Yamnaya", "French"]


_R_QPGRAPH = """
suppressMessages(library(admixtools))
# admixtools2 2.0.10 qpgraph() with igraph >= 2, as in tests/golden/admixtools2/generate_reference.R
local({
  imports <- parent.env(asNamespace("admixtools"))
  shim <- function(graph, vp, ...) igraph::get_edge_ids(graph, as.vector(vp), ...)
  unlockBinding("get.edge.ids", imports)
  assign("get.edge.ids", shim, envir = imports)
  lockBinding("get.edge.ids", imports)
})
a <- commandArgs(TRUE)
edges <-as.matrix(read.table(a[2], col.names = c("from", "to")))
f2 <- f2_from_geno(a[1], pops = strsplit(a[3], ",")[[1]], verbose = FALSE)
r <- qpgraph(f2, edges, numstart = 10, seed = as.integer(a[4]), verbose = FALSE)
jsonlite::write_json(list(score = r$score), file.path(a[5], "qpgraph.json"), digits = NA, auto_unbox = TRUE)
"""


def _graph_file(o: Path) -> Path:
    path = o / "graph.txt"
    path.write_text("".join(f"{a}\t{b}\n" for a, b in GRAPH), encoding="utf-8")
    return path


def _score(name: str) -> Callable[[Path], float]:
    return lambda d: json.loads((d / name).read_text(encoding="utf-8"))["score"]


def _popmix_qpgraph(i: Path, o: Path, t: int, seed: int = 1) -> list[str]:
    return ["sh", "-c", f"popmix --threads {t} qpgraph -i {i} --graph {_graph_file(o)} "
                        f"--numstart 10 --seed {seed} --format json > {o / 'qpgraph.json'}"]


def _ref_qpgraph(i: Path, o: Path, seed: int = 1) -> list[str]:
    return _r(_R_QPGRAPH, o, str(i), str(_graph_file(o)), ",".join(GRAPH_POPS), str(seed))


def _agree_qpgraph(ours: Path, theirs: Path) -> dict:
    got, want = _score("qpgraph.json")(ours), _score("qpgraph.json")(theirs)
    rel = abs(got - want) / abs(want)
    # the best of 10 random starts, by seed: each tool's run-to-run spread
    mine = [got, *_seed_runs(theirs, "popmix", lambda d, s: _popmix_qpgraph(subset(), d, 4, s),
                             _score("qpgraph.json")).values()]
    ref = [want, *_seed_runs(theirs, "admixtools2", lambda d, s: _ref_qpgraph(subset(), d, s),
                             _score("qpgraph.json")).values()]
    # an optimized score: R's L-BFGS-B stops at factr * epsilon
    return {"score": want, "score_rel": rel, "seeds": 1 + len(SPREAD_SEEDS),
            "popmix_spread": max(mine) - min(mine), "reference_spread": max(ref) - min(ref),
            "pass": rel <= 1e-6}


QPGRAPH_CASE = Case(
    name="qpgraph", analysis="qpGraph, 6 populations, 1 admixture, 10 starts",
    reference="admixtools2 2.0.10 qpgraph",
    popmix=_popmix_qpgraph,
    ref=_ref_qpgraph,
    agree=_agree_qpgraph,
)


FIND_GRAPHS_ADMIX = 1


# admixtools2's evolutionary search: generations as in its defaults, fewer of them
FIND_GRAPHS_GENERATIONS = 20


_R_FIND_GRAPHS = """
suppressMessages(library(admixtools))
future::plan("sequential")
local({
  imports <- parent.env(asNamespace("admixtools"))
  shim <- function(graph, vp, ...) igraph::get_edge_ids(graph, as.vector(vp), ...)
  unlockBinding("get.edge.ids", imports)
  assign("get.edge.ids", shim, envir = imports)
  lockBinding("get.edge.ids", imports)
})
a <- commandArgs(TRUE)
pops <- strsplit(a[2], ",")[[1]]
f2 <- f2_from_geno(a[1], pops = pops, verbose = FALSE)
set.seed(as.integer(a[5]))
r <- find_graphs(f2, numadmix = as.integer(a[3]), outpop = pops[1], stop_gen = as.integer(a[4]), verbose = FALSE)
jsonlite::write_json(list(score = min(r$score), graphs = nrow(r)), file.path(a[6], "find_graphs.json"),
                     digits = NA, auto_unbox = TRUE)
"""


def _ref_find_graphs(i: Path, o: Path, seed: int = 1) -> list[str]:
    return _r(_R_FIND_GRAPHS, o, str(i), ",".join(GRAPH_POPS), str(FIND_GRAPHS_ADMIX),
              str(FIND_GRAPHS_GENERATIONS), str(seed))


def _agree_find_graphs(ours: Path, theirs: Path) -> dict:
    got = json.loads((ours / "find_graphs.json").read_text(encoding="utf-8"))
    want = _score("find_graphs.json")(theirs)
    best = got["graphs"][0]["score"]
    # the evolutionary search's best by seed; popmix's exhaustive search has
    # no seed
    ref = [want, *_seed_runs(theirs, "admixtools2", lambda d, s: _ref_find_graphs(subset(), d, s),
                             _score("find_graphs.json")).values()]
    # an exhaustive search: its best is at least as good as any graph's,
    # each fit to L-BFGS-B's convergence
    return {"graphs_fitted": got["n_evaluated"], "score": best, "reference_score": want,
            "score_rel": (best - want) / abs(want), "seeds": 1 + len(SPREAD_SEEDS),
            "reference_spread": max(ref) - min(ref), "pass": best <= want * (1 + 1e-6)}


FIND_GRAPHS_CASE = Case(
    name="find_graphs", analysis=f"find_graphs, {len(GRAPH_POPS)} populations, {FIND_GRAPHS_ADMIX} admixture, "
                                 "Chimp as outgroup: popmix exhaustive, admixtools2 evolutionary "
                                 f"({FIND_GRAPHS_GENERATIONS} generations, seed 1)",
    reference="admixtools2 2.0.10 find_graphs",
    popmix=lambda i, o, t: ["sh", "-c", f"popmix --threads {t} find-graphs -i {i} --populations "
                                        f"{','.join(GRAPH_POPS)} --nadmix {FIND_GRAPHS_ADMIX} --outgroup "
                                        f"{GRAPH_POPS[0]} --method exhaustive --top 1 --format json "
                                        f"> {o / 'find_graphs.json'}"],
    ref=_ref_find_graphs,
    agree=_agree_find_graphs,
)
