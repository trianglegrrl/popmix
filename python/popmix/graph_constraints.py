"""Constraints on admixture graphs, as admixtools2's ``summarize_numadmix``,
``satisfies_numadmix``, ``summarize_zerof4``, ``satisfies_zerof4``,
``satisfies_nonzerof4``, ``satisfies_eventorder`` and
``satisfies_constraints``.

A graph is a list of (from, to) edges, as for ``qpgraph``. ``find_graphs``
uses these to restrict its search.
"""

from __future__ import annotations

import itertools
import random
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

Edges = Sequence[tuple[str, str]]


class _Graph:
    def __init__(self, edges: Edges):
        self.edges = [tuple(e) for e in edges]
        self.parents: dict[str, list[str]] = {}
        self.children: dict[str, list[str]] = {}
        nodes: list[str] = []
        for a, b in self.edges:
            for v in (a, b):
                if v not in self.parents:
                    self.parents[v], self.children[v] = [], []
                    nodes.append(v)
            self.parents[b].append(a)
            self.children[a].append(b)
        roots = [v for v in nodes if not self.parents[v]]
        if len(roots) != 1:
            raise ValueError(f"the graph needs exactly one root, found {roots}")
        self.root = roots[0]
        self.leaves = [v for v in nodes if not self.children[v]]
        self.admix = [v for v in nodes if len(self.parents[v]) > 1]

    def root_paths(self, leaf: str) -> list[list[str]]:
        """Every root-to-leaf path, as vertex lists."""
        out: list[list[str]] = []

        def walk(v: str, path: list[str]):
            if v == self.root:
                out.append(path[::-1])
                return
            for p in self.parents[v]:
                walk(p, path + [p])

        walk(leaf, [leaf])
        return out

    def up_paths(self, v: str) -> list[list[str]]:
        """Every path from `v` up to the root (igraph all_simple_paths mode 'in')."""
        return [p[::-1] for p in self.root_paths(v)]


def summarize_numadmix(edges: Edges) -> dict[str, int]:
    """The largest number of admixture events on a path from the root to each leaf."""
    g = _Graph(edges)
    adm = set(g.admix)
    return {leaf: max(sum(v in adm for v in p) for p in g.root_paths(leaf)) for leaf in g.leaves}


def _bound(x: Any) -> float | None:
    if x is None or (isinstance(x, float) and x != x):
        return None
    return x


def satisfies_numadmix(edges: Edges, admix_constraints: Iterable[Mapping[str, Any]] | None) -> bool:
    """Whether each population's number of admixture events is within its
    constraint.

    ``admix_constraints`` lists ``{"pop", "min", "max"}``; a missing (None or
    NaN) bound is no bound, and populations not in the graph are ignored.
    """
    if admix_constraints is None:
        return True
    n = summarize_numadmix(edges)
    for c in admix_constraints:
        if not {"pop", "min", "max"} <= set(c):
            raise ValueError("'admix_constraints' should have columns 'pop', 'min', 'max'!")
        k = n.get(c["pop"])
        lo, hi = _bound(c["min"]), _bound(c["max"])
        if k is not None and ((lo is not None and k < lo) or (hi is not None and k > hi)):
            return False
    return True


def _f2_path_pairs(g: _Graph, rng: random.Random) -> dict[tuple[str, str], float]:
    """admixtools2's expected f2 for random edge lengths and admixture weights:
    over every pair of root-to-leaf paths, the product of their admixture
    weights times the length of the edges on exactly one of them."""
    weight = {v: rng.random() for v in g.admix}
    length = {e: rng.random() for e in g.edges}
    paths = {}
    for leaf in g.leaves:
        ps = []
        for p in g.root_paths(leaf):
            w = 1.0
            for a, b in zip(p, p[1:]):  # noqa: B905
                if b in weight:
                    w *= weight[b] if g.parents[b][0] == a else 1 - weight[b]
            ps.append((w, set(zip(p, p[1:]))))  # noqa: B905
        paths[leaf] = ps
    f2 = {}
    for a, b in itertools.permutations(g.leaves, 2):
        f2[a, b] = sum(wa * wb * sum(length[e] for e in ea ^ eb) for wa, ea in paths[a] for wb, eb in paths[b])
    return f2


def summarize_zerof4(edges: Edges, eps: float = 1e-7, seed: int = 1) -> set[tuple[str, str, str, str]]:
    """The f4-statistics the graph's topology makes zero, as admixtools2
    ``summarize_zerof4``: f4(pop1, pop2; pop3, pop4) with pop1 the first of
    the four and pop1 < pop2, pop3 < pop4, zero for random edge lengths and
    admixture weights."""
    g = _Graph(edges)
    f2 = _f2_path_pairs(g, random.Random(seed))
    pops = g.leaves
    zero = set()
    for p1, p2, p3, p4 in itertools.product(pops, repeat=4):
        if p1 < p2 and p1 < p3 and p1 < p4 and p3 < p4 and p2 != p3 and p2 != p4:
            f4 = (f2[p1, p4] + f2[p2, p3] - f2[p1, p3] - f2[p2, p4]) / 2
            if -eps <= f4 <= eps:
                zero.add((p1, p2, p3, p4))
    return zero


def _normalize_f4(quads: Iterable[Sequence[str]]) -> set[tuple[str, str, str, str]]:
    out = set()
    for q in quads:
        a, b, c, d = q[:4]
        p1, p2, p3, p4 = min(a, b), max(a, b), min(c, d), max(c, d)
        out.add((p1, p2, p3, p4) if p1 < p3 else (p3, p4, p1, p2))
    return out


def satisfies_zerof4(edges: Edges, zero_f4: Iterable[Sequence[str]] | None) -> bool:
    """Whether every f4-statistic in ``zero_f4`` (4 populations each) is zero
    in the graph."""
    want = _normalize_f4(zero_f4 or [])
    return not want or want <= summarize_zerof4(edges)


def satisfies_nonzerof4(edges: Edges, nonzero_f4: Iterable[Sequence[str]] | None) -> bool:
    """Whether no f4-statistic in ``nonzero_f4`` is zero in the graph."""
    want = _normalize_f4(nonzero_f4 or [])
    return not want or not (want & summarize_zerof4(edges))


def _intersect(a: list[str], b: list[str]) -> list[str]:
    s = set(b)
    return [v for v in dict.fromkeys(a) if v in s]


def _setdiff(a: list[str], b: list[str]) -> list[str]:
    s = set(b)
    return [v for v in dict.fromkeys(a) if v not in s]


def _one_event(g: _Graph, earlier1: str, earlier2: str | None, later1: str, later2: str | None,
               type_: int) -> bool | None:
    if earlier2 is None:
        (earlier2,) = g.parents[earlier1]
        if earlier2 == g.root:
            return True
    if later2 is None:
        (later2,) = g.parents[later1]
        if later2 == g.root:
            return False
    pe1, pe2, pl1, pl2 = (g.up_paths(v) for v in (earlier1, earlier2, later1, later2))
    left = [_intersect(a, b) for a in pe1 for b in pe2]
    right = [_intersect(a, b) for a in pl1 for b in pl2]
    if type_ != 1:
        return all(is1[0] in is2 for is1 in left for is2 in right)
    d1 = [[len(_setdiff(is1, is2)) for is1 in left] for is2 in right]  # per j, over i
    d2 = [[len(_setdiff(is2, is1)) for is2 in right] for is1 in left]  # per i, over j
    m10 = min(max(row) for row in d1) == 0
    m20 = min(max(row) for row in d2) == 0
    if m10 and not m20:
        return True
    if m20 and not m10:
        return False
    return None


def satisfies_eventorder(edges: Edges, eventorder: Iterable[Mapping[str, Any]] | None, strict: bool = True) -> bool:
    """Whether the graph puts each split of ``earlier1`` and ``earlier2``
    before that of ``later1`` and ``later2`` (admixtools2
    ``satisfies_eventorder``).

    Each constraint is ``{"earlier1", "earlier2", "later1", "later2"}`` and
    optionally ``"type"``: 1 (default), some lineage puts the earlier split
    first; 2, every lineage does. A missing ``earlier2`` or ``later2`` (None)
    stands for the parent of ``earlier1`` or ``later1``. A constraint the
    graph does not decide fails with ``strict`` and is ignored without.
    """
    if eventorder is None:
        return True
    g = _Graph(edges)
    status = [_one_event(g, c["earlier1"], _bound(c.get("earlier2")), c["later1"], _bound(c.get("later2")),
                         int(c.get("type", 1))) for c in eventorder]
    if strict:
        return all(s is True for s in status)
    return all(s for s in status if s is not None)


def satisfies_constraints(edges: Edges, nonzero_f4=None, admix_constraints=None, event_order=None) -> bool:
    """All of ``satisfies_nonzerof4``, ``satisfies_numadmix`` and ``satisfies_eventorder``."""
    return (satisfies_nonzerof4(edges, nonzero_f4) and satisfies_numadmix(edges, admix_constraints)
            and satisfies_eventorder(edges, event_order))
