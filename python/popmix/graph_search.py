"""admixtools2's evolutionary graph search (``find_graphs``).

``find_graphs(method="evolutionary")`` searches as admixtools2's
``find_graphs`` does. It starts from a random graph with ``nadmix``
admixture nodes, then in each generation takes the best graph found so far
and makes up to ``numgraphs`` new ones from it:

- for its most negative edges in an unconstrained fit, a rearrangement that
  moves them (admixtools2 ``rearrange_negadmix3`` or
  ``replace_admix_with_random`` for admixture edges, ``rearrange_negdrift``
  for drift edges), then swaps of pairs of negative edges
  (``swap_negdrift``);
- the rest by random mutations: subtree prune and regraft of leaves or of
  any subtree, swapping two leaves, moving an admixture edge, flipping one,
  moving the root, and two of these at once (``spr_leaves``, ``spr_all``,
  ``swap_leaves``, ``move_admixedge_once``, ``flipadmix_random``,
  ``place_root_random``, ``mutate_n``).

Graphs already seen (by ``graph_key``), graphs without the outgroup as a
child of the root, and graphs that break a constraint once some graph has
met it are dropped; the rest are fitted with ``qpgraph`` from one random
start. Every ``plusminus_generations`` generations without improvement it
instead adds an admixture edge to the best graph (the best of up to ten
times ``numgraphs`` random ones) and, beyond ``max_admix``, removes the one
that fits worst; the result is kept only if it is at least 5% better. It
stops after ``stop_gen`` generations, ``stop_gen2`` without improvement,
at ``stop_score`` or after ``stop_sec`` seconds.

The operations are admixtools2's, applied to edge lists; random choices
come from popmix's own generator (``seed``), so the path differs from
admixtools2's for the same data.
"""

from __future__ import annotations

import itertools
import random
import time
from collections import deque
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

from popmix import _core
from popmix.dstats import qpdstat
from popmix.graph_constraints import satisfies_eventorder, satisfies_nonzerof4, satisfies_numadmix

Edge = tuple[str, str]


@dataclass
class EvolvedGraph:
    """One graph the evolutionary search fitted."""

    edges: list[Edge]
    """The graph, the path from the root to the f3 base population first."""
    score: float
    generation: int
    mutation: str


# ---------------------------------------------------------------------------
# edge-list graphs


def _nodes(edges: Sequence[Edge]) -> list[str]:
    return list(dict.fromkeys(v for e in edges for v in e))


def _parents(edges: Sequence[Edge], v: str) -> list[str]:
    return [a for a, b in edges if b == v]


def _children(edges: Sequence[Edge], v: str) -> list[str]:
    return [b for a, b in edges if a == v]


def _root(edges: Sequence[Edge]) -> str:
    targets = {b for _, b in edges}
    roots = [v for v in _nodes(edges) if v not in targets]
    if len(roots) != 1:
        raise ValueError("not one root")
    return roots[0]


def _leaves(edges: Sequence[Edge]) -> list[str]:
    sources = {a for a, _ in edges}
    return [v for v in _nodes(edges) if v not in sources]


def _admix_nodes(edges: Sequence[Edge]) -> list[str]:
    return [v for v in _nodes(edges) if len(_parents(edges, v)) == 2]


def _numadmix(edges: Sequence[Edge]) -> int:
    return len(_admix_nodes(edges))


def _descendants(edges: Sequence[Edge], v: str) -> set[str]:
    out, todo = {v}, [v]
    while todo:
        for c in _children(edges, todo.pop()):
            if c not in out:
                out.add(c)
                todo.append(c)
    return out


def _outpop(edges: Sequence[Edge]) -> str | None:
    """The leaf directly below the root, if there is exactly one (admixtools2 ``get_outpop``)."""
    leaves = set(_leaves(edges))
    out = [c for c in _children(edges, _root(edges)) if c in leaves]
    return out[0] if len(out) == 1 else None


def _is_valid(edges: Sequence[Edge]) -> bool:
    """admixtools2 ``is_valid``: simple, connected, acyclic, one root with two
    children, in- and out-degree at most 2, leaves with one parent."""
    if len(set(edges)) != len(edges) or any(a == b for a, b in edges):
        return False
    nodes = _nodes(edges)
    indeg = dict.fromkeys(nodes, 0)
    outdeg = dict.fromkeys(nodes, 0)
    for a, b in edges:
        outdeg[a] += 1
        indeg[b] += 1
    roots = [v for v in nodes if indeg[v] == 0]
    if len(roots) != 1 or outdeg[roots[0]] < 2:
        return False
    if any(indeg[v] > 2 or outdeg[v] > 2 or (outdeg[v] == 0 and indeg[v] != 1) for v in nodes):
        return False
    # acyclic and every node reachable from the root (Kahn)
    left = dict(indeg)
    todo = deque(roots)
    seen = 0
    while todo:
        v = todo.popleft()
        seen += 1
        for c in _children(edges, v):
            left[c] -= 1
            if left[c] == 0:
                todo.append(c)
    return seen == len(nodes)


class _Names:
    def __init__(self, edges: Sequence[Edge]):
        self.taken = set(_nodes(edges))
        self.n = 0

    def new(self) -> str:
        while True:
            self.n += 1
            name = f"n{self.n}"
            if name not in self.taken:
                self.taken.add(name)
                return name


def _simplify(edges: Sequence[Edge], names: _Names) -> list[Edge]:
    """admixtools2 ``simplify_graph``: no repeated edges, no node with one
    parent and one child, and nodes with two parents and two children split
    into an admixture node and a split node (``X_to_H``)."""
    edges = list(dict.fromkeys(edges))
    while True:
        changed = False
        for v in _nodes(edges):
            ps, cs = _parents(edges, v), _children(edges, v)
            if len(ps) == 1 and len(cs) == 1:
                edges = [e for e in edges if v not in e] + [(ps[0], cs[0])]
                edges = list(dict.fromkeys(edges))
                changed = True
                break
            if len(ps) > 1 and len(cs) > 1:
                h = names.new()
                edges = [e for e in edges if e[1] != v] + [(p, h) for p in ps] + [(h, v)]
                changed = True
                break
        if not changed:
            return edges


def _prune(edges: list[Edge], leaves: set[str], names: _Names) -> list[Edge]:
    """Remove internal nodes that became leaves, and a root with one child."""
    edges = _simplify(edges, names)
    while True:
        extra = [v for v in _leaves(edges) if v not in leaves]
        if not extra:
            break
        edges = _simplify([e for e in edges if e[1] not in extra], names)
    root = _root(edges)
    if len(_children(edges, root)) == 1:
        edges = [e for e in edges if e[0] != root]
    return edges


def _ok(edges: Sequence[Edge] | None, leaves: Sequence[str]) -> bool:
    return edges is not None and _is_valid(edges) and sorted(_leaves(edges)) == sorted(leaves)


def base_first(edges: Sequence[Edge], base: str) -> list[Edge]:
    """The edges with a path from the root to ``base`` first, so that
    ``qpgraph`` takes ``base`` as its f3 base population."""
    path, v = [], base
    while _parents(edges, v):
        p = _parents(edges, v)[0]
        path.append((p, v))
        v = p
    path.reverse()
    return path + [e for e in edges if e not in set(path)]


# ---------------------------------------------------------------------------
# admixtools2's graph operations


def _new_edges(edges: Sequence[Edge], fix_outgroup: bool = True) -> list[tuple[str, str, str, str]]:
    """admixtools2 ``find_newedges``: pairs of edges an admixture edge can
    join, the second not above the first."""
    cand = list(edges)
    if fix_outgroup:
        out = _outpop(edges)
        if out is not None:
            cand = [e for e in cand if e != (_root(edges), out)]
    below = {v: _descendants(edges, v) for v in _nodes(edges)}
    return [(sf, st, df, dt) for (sf, st), (df, dt) in itertools.product(cand, cand)
            if sf != df and st != dt and st != df and dt != sf and sf not in below[dt]]


def _insert_admix(edges: Sequence[Edge], names: _Names, rng: random.Random, new: tuple | None = None,
                  fix_outgroup: bool = True) -> list[Edge] | None:
    """admixtools2 ``insert_admix``: an admixture edge from a new node on
    sf -> st to a new node on df -> dt (a random such pair if not given or
    not possible)."""
    leaves = _leaves(edges)
    if new is not None:
        sf, st, df, dt = new
        if (sf, st) in edges and (df, dt) in edges:
            s, d = names.new(), names.new()
            g = [e for e in edges if e not in {(sf, st), (df, dt)}] + [(sf, s), (s, st), (df, d), (d, dt), (s, d)]
            if _ok(g, leaves):
                return g
    cands = _new_edges(edges, fix_outgroup)
    rng.shuffle(cands)
    for c in cands[:20]:
        g = _insert_admix(edges, names, rng, c, fix_outgroup) if c != new else None
        if g is not None:
            return g
    return None


def _insert_admix_above(edges: Sequence[Edge], names: _Names, rng: random.Random, source_to: str,
                        dest_to: str) -> list[Edge] | None:
    """An admixture edge from above ``source_to`` to above ``dest_to``, or a random one."""
    sp, dp = _parents(edges, source_to), _parents(edges, dest_to)
    new = (rng.choice(sp), source_to, rng.choice(dp), dest_to) if sp and dp else None
    return _insert_admix(edges, names, rng, new)


def _delete_admix(edges: Sequence[Edge], names: _Names, frm: str, to: str) -> list[Edge] | None:
    """admixtools2 ``delete_admix``: the graph without the admixture edge frm -> to."""
    leaves = _leaves(edges)
    g = _prune([e for e in edges if e != (frm, to)], set(leaves), names)
    return g if _ok(g, leaves) else None


def _split_graph(edges: Sequence[Edge], names: _Names, rng: random.Random):
    """admixtools2 ``split_graph``: a tree, removing one incoming edge of each
    admixture node, and the (source_to, dest_to) pairs that restore them."""
    tree = list(edges)
    restore = []
    leaves = set(_leaves(edges))
    while True:
        adm = _admix_nodes(tree)
        if not adm:
            return tree, restore
        a = adm[0]
        p = rng.choice(_parents(tree, a))
        sibling = [c for c in _children(tree, p) if c != a]
        child = _children(tree, a)
        restore.append((sibling[0] if sibling else p, child[0] if child else a))
        tree = _prune([e for e in tree if e != (p, a)], leaves, names)


def _spr_tree(edges: Sequence[Edge], names: _Names, rng: random.Random, only_leaves: bool,
              fix_outgroup: bool) -> list[Edge] | None:
    """admixtools2 ``subtree_prune_and_regraft`` on a tree."""
    root = _root(edges)
    firstgen = set(_children(edges, root)) if fix_outgroup else set()
    leaves = set(_leaves(edges))
    nodes = _nodes(edges)
    for _ in range(100):
        cand = [v for v in nodes if v != root and v not in firstgen and (not only_leaves or v in leaves)]
        if not cand:
            return None
        cut = rng.choice(cand)
        sub = _descendants(edges, cut)
        (cp,) = _parents(edges, cut)
        gp = _parents(edges, cp)
        sib = [c for c in _children(edges, cp) if c != cut]
        hosts = [v for v in nodes if v not in {root, cp} | sub | set(sib) | firstgen]
        if hosts and sib and gp:
            break
    else:
        return None
    host = rng.choice(hosts)
    (hp,) = _parents(edges, host)
    n = names.new()
    g = [e for e in edges if cp not in e and e != (hp, host)]
    g += [(gp[0], sib[0]), (hp if hp != cp else gp[0], n), (n, host), (n, cut)]
    return g


def spr(edges: Sequence[Edge], names: _Names, rng: random.Random, only_leaves: bool = False,
        fix_outgroup: bool = True) -> list[Edge] | None:
    """admixtools2 ``spr_leaves`` / ``spr_all``: prune a leaf (or any subtree)
    of the graph's tree and regraft it elsewhere, then restore the admixture
    edges where they were, or at random where they no longer fit."""
    leaves = _leaves(edges)
    tree, restore = _split_graph(edges, names, rng)
    g = _spr_tree(tree, names, rng, only_leaves, fix_outgroup)
    for source_to, dest_to in restore:
        if g is None:
            return None
        g = _insert_admix_above(g, names, rng, source_to, dest_to)
    return g if _ok(g, leaves) else None


def swap_leaves(edges: Sequence[Edge], names: _Names, rng: random.Random, fix_outgroup: bool = True):
    """admixtools2 ``swap_leaves``: exchange two leaves (not the first, the outgroup)."""
    leaves = _leaves(edges)
    pool = leaves[1:] if fix_outgroup else leaves
    if len(pool) < 2:
        return None
    a, b = rng.sample(pool, 2)
    swap = {a: b, b: a}
    return [(swap.get(x, x), swap.get(y, y)) for x, y in edges]


def _bfs(edges: Sequence[Edge], start: str) -> list[str]:
    seen, order, todo = {start}, [start], deque([start])
    while todo:
        v = todo.popleft()
        for w in _parents(edges, v) + _children(edges, v):
            if w not in seen:
                seen.add(w)
                order.append(w)
                todo.append(w)
    return order


def move_admixedge_once(edges: Sequence[Edge], names: _Names, rng: random.Random, fix_outgroup: bool = True):
    """admixtools2 ``move_admixedge_once``: move one end of an admixture edge
    to the nearest place it fits."""
    leaves = _leaves(edges)
    root = _root(edges)
    firstgen = set(_children(edges, root)) if fix_outgroup else set()
    admix = _admix_nodes(edges)
    rng.shuffle(admix)
    for a in admix:
        parents = _parents(edges, a)
        rng.shuffle(parents)
        below_a = _descendants(edges, a)
        near = set(_children(edges, parents[0])) | set(_children(edges, parents[1]))
        for p in parents:
            sib = [c for c in _children(edges, p) if c != a]
            gp = _parents(edges, p)
            if not sib or not gp:
                continue
            for n in _bfs(edges, p):
                if (n == root or n in firstgen or n in parents or n in gp or n in admix or n in near
                        or n in below_a):
                    continue
                ngp = _parents(edges, n)
                g = [e for e in edges if e not in {(ngp[0], n), (gp[0], p), (p, sib[0])}]
                g += [(gp[0], sib[0]), (ngp[0], p), (p, n)]
                g = _simplify(g, names)
                return g if _ok(g, leaves) else None
    return None


def flipadmix_random(edges: Sequence[Edge], names: _Names, rng: random.Random, fix_outgroup: bool = True):
    """admixtools2 ``flipadmix_random``: reverse an admixture edge (and the
    edge above it, if that gives a valid graph)."""
    leaves = _leaves(edges)
    adm = [(p, a) for a in _admix_nodes(edges) for p in _parents(edges, a)]
    rng.shuffle(adm)
    for frm, to in adm:
        g = [(to, frm) if e == (frm, to) else e for e in edges]
        above = [k for k, e in enumerate(g) if e[1] == frm]
        g2 = list(g)
        if above:
            x, y = g2[above[0]]
            g2[above[0]] = (y, x)
        for cand in (g2, g):
            if _ok(cand, leaves):
                return cand
    return None


def place_root_random(edges: Sequence[Edge], names: _Names, rng: random.Random, fix_outgroup: bool = True):
    """admixtools2 ``place_root_random``: put the root on a random edge (below
    the outgroup's split, which stays at the top)."""
    leaves = _leaves(edges)
    root = _root(edges)
    out = _outpop(edges) if fix_outgroup else None
    g = list(edges)
    if out is not None:
        g = [e for e in g if root not in e]
        root = _root(g)
    inner = [e for e in g if root not in e and e[1] not in set(_leaves(g))]
    if not inner:
        return None
    frm, to = rng.choice(inner)
    # remove the root, join its children, root the result on frm -> to
    kids = _children(g, root)
    und = [e for e in g if root not in e and e != (frm, to)]
    if len(kids) == 2:
        und.append((kids[0], kids[1]))
    newroot = names.new()
    und += [(newroot, frm), (newroot, to)]
    nbrs: dict[str, list[str]] = {}
    for a, b in und:
        nbrs.setdefault(a, []).append(b)
        nbrs.setdefault(b, []).append(a)
    dist, todo = {newroot: 0}, deque([newroot])
    while todo:
        v = todo.popleft()
        for w in nbrs[v]:
            if w not in dist:
                dist[w] = dist[v] + 1
                todo.append(w)
    oriented = [(a, b) if dist[a] <= dist[b] else (b, a) for a, b in und]
    keep = set(leaves) - ({out} if out else set())
    try:
        res = _prune(oriented, keep, names)
    except ValueError:
        return None
    if out is not None:
        top = names.new()
        res = [(top, out), (top, _root(res))] + res
    return res if _ok(res, leaves) else None


def mutate_n(edges, names, rng, fix_outgroup=True, n: int = 2):
    """admixtools2 ``mutate_n``: two random mutations in a row."""
    funs = [lambda *a, **k: spr(*a, **k), lambda *a, **k: spr(*a, only_leaves=True, **k), swap_leaves,
            move_admixedge_once, flipadmix_random, mutate_n]
    g = list(edges)
    for f in rng.choices(funs, k=n):
        g = f(g, names, rng, fix_outgroup=fix_outgroup)
        if g is None:
            return None
    return g


MUTATIONS: dict[str, Callable] = {
    "spr_leaves": lambda g, nm, r, fix_outgroup=True: spr(g, nm, r, only_leaves=True, fix_outgroup=fix_outgroup),
    "spr_all": lambda g, nm, r, fix_outgroup=True: spr(g, nm, r, only_leaves=False, fix_outgroup=fix_outgroup),
    "swap_leaves": swap_leaves,
    "move_admixedge_once": move_admixedge_once,
    "flipadmix_random": flipadmix_random,
    "place_root_random": place_root_random,
    "mutate_n": mutate_n,
}


def rearrange_negdrift(edges, names, rng, frm: str, to: str):
    """admixtools2 ``rearrange_negdrift``: swap a negative drift edge's lower
    end with its sibling."""
    sib = [c for c in _children(edges, frm) if c != to]
    kids = _children(edges, to)
    if not kids or not sib:
        return list(edges)
    c = rng.choice(kids)
    g = [e for e in edges if e not in {(frm, sib[0]), (to, c)}] + [(to, sib[0]), (frm, c)]
    if _ok(g, _leaves(edges)) and _outpop(g) == _outpop(edges):
        return g
    return None


def rearrange_negadmix3(edges, names, rng, frm: str, to: str):
    """admixtools2 ``rearrange_negadmix3``: move an admixed node with a
    negative admixture edge above its other parent."""
    leaves = _leaves(edges)
    pos = [p for p in _parents(edges, to) if p != frm]
    if len(pos) != 1:
        return None
    gp = _parents(edges, pos[0])
    if len(gp) != 1:
        return None
    g = [e for e in edges if e not in {(gp[0], pos[0]), (pos[0], to), (frm, to)}] + [(gp[0], to), (to, pos[0])]
    if not _is_valid(g):
        return None
    g2 = _simplify(g + [(frm, pos[0])], names)
    if _is_valid(g2):
        g = g2
    else:
        g = _insert_admix(g, names, rng)
        if g is None:
            return None
        g = _simplify(g, names)
    if not _ok(g, leaves):
        return None
    for _ in range(_numadmix(edges) - _numadmix(g)):
        g = _insert_admix(g, names, rng)
        if g is None:
            return None
    return g if _ok(g, leaves) else None


def replace_admix_with_random(edges, names, rng, frm: str, to: str):
    """admixtools2 ``replace_admix_with_random``: remove an admixture edge and add a random one."""
    g = _delete_admix(edges, names, frm, to)
    return None if g is None else _insert_admix(g, names, rng)


def swap_negdrift(edges, names, rng, negative: Sequence[Edge]) -> list[list[Edge]]:
    """admixtools2 ``swap_negdrift``: for each pair of negative edges, the
    graph with their lower ends exchanged."""
    out = []
    leaves = _leaves(edges)
    for (f1, t1), (f2, t2) in itertools.product(negative, negative):
        if f1 != f2 and t1 != t2 and f1 != t2 and t1 != f2:
            g = [e for e in edges if e not in {(f1, t1), (f2, t2)}] + [(f1, t2), (f2, t1)]
            if _ok(g, leaves):
                out.append(g)
    return out


def random_admixturegraph(pops: Sequence[str], nadmix: int, rng: random.Random, outgroup: str | None = None,
                          ntry: int = 100) -> list[Edge]:
    """A random graph with these leaves and admixture nodes (admixtools2
    ``random_admixturegraph``), with ``outgroup``, if given, below the root."""
    names = _Names([("R", p) for p in pops])
    rest = [p for p in pops if p != outgroup]
    for _ in range(ntry):
        order = rest[:]
        rng.shuffle(order)
        top = names.new()
        edges = [(top, order[0]), (top, order[1])]
        for leaf in order[2:]:
            a, b = rng.choice(edges)
            n = names.new()
            edges = [e for e in edges if e != (a, b)] + [(a, n), (n, b), (n, leaf)]
        if outgroup is not None:
            r = names.new()
            edges = [(r, outgroup), (r, top)] + edges
        g = edges
        for _ in range(nadmix):
            g = _insert_admix(g, names, rng, fix_outgroup=outgroup is not None)
            if g is None:
                break
        if g is not None and _numadmix(g) == nadmix:
            return g
    raise ValueError(f"could not make a graph with {nadmix} admixture nodes")


# ---------------------------------------------------------------------------
# the search


def evolve_graphs(
    f2_blocks: _core.PyF2Blocks,
    populations: Sequence[str],
    nadmix: int = 0,
    outgroup: str | None = None,
    *,
    initgraph: Sequence[Edge] | None = None,
    numgraphs: int = 10,
    stop_gen: int = 100,
    stop_gen2: int = 15,
    stop_score: float = 0.0,
    stop_sec: float | None = None,
    plusminus_generations: int = 5,
    max_admix: int | None = None,
    admix_constraints: Sequence[dict] | None = None,
    event_constraints: Sequence[dict] | None = None,
    reject_f4z: float = 0.0,
    seed: int = 1,
    numstart: int = 1,
    diag: float = 1e-4,
    diag_f3: float = 1e-5,
) -> tuple[list[EvolvedGraph], int]:
    """See ``find_graphs(method="evolutionary")``."""
    rng = random.Random(seed)
    pops = list(populations)
    base = pops[0]
    fix_outgroup = outgroup is not None
    graph = [tuple(e) for e in initgraph] if initgraph is not None else \
        random_admixturegraph(pops, nadmix, rng, outgroup)
    if sorted(_leaves(graph)) != sorted(pops):
        raise ValueError("the initial graph's leaves must be the populations")
    nadmix = _numadmix(graph)
    max_admix = nadmix if max_admix is None else max_admix
    names = _Names(graph)

    nonzero_f4 = None
    if reject_f4z > 0:
        t = qpdstat(f2_blocks, pops, pops, pops, pops)
        nonzero_f4 = [tuple(r) for r in t.loc[t.z.abs() > reject_f4z, ["pop1", "pop2", "pop3", "pop4"]].itertuples(
            index=False)]

    def fit(g: Sequence[Edge], **kw: Any) -> Any:
        return _core.qpgraph(f2_blocks, base_first(g, base), numstart=numstart, seed=rng.randrange(1, 2**31 - 1),
                             diag=diag, diag_f3=diag_f3, **kw)

    def key(g: Sequence[Edge]) -> tuple:
        return tuple(_core.graph_key(list(g)))

    models: dict[tuple, EvolvedGraph] = {}

    def add(g: list[Edge], gen: int, mut: str) -> EvolvedGraph:
        m = EvolvedGraph(base_first(g, base), fit(g).score, gen, mut)
        models[key(g)] = m
        return m

    best = add(_simplify(graph, names), 0, "none").score
    gimp = 0
    active = {"nzf4": None, "admix": None, "event": None}

    def allowed(g: list[Edge]) -> bool:
        return (satisfies_nonzerof4(g, active["nzf4"]) and satisfies_numadmix(g, active["admix"])
                and satisfies_eventorder(g, active["event"]))

    stop_at = None if stop_sec is None else time.monotonic() + stop_sec
    for gen in range(1, stop_gen + 1):
        if best <= stop_score or gimp > stop_gen2 or (stop_at is not None and time.monotonic() > stop_at):
            break
        sel = min(models.values(), key=lambda m: m.score)
        graph = sel.edges
        if active["nzf4"] is None and nonzero_f4 and satisfies_nonzerof4(graph, nonzero_f4):
            active["nzf4"] = nonzero_f4
        if active["admix"] is None and admix_constraints and satisfies_numadmix(graph, admix_constraints):
            active["admix"] = admix_constraints
        if active["event"] is None and event_constraints and satisfies_eventorder(graph, event_constraints):
            active["event"] = event_constraints

        if gimp > 0 and gimp % plusminus_generations == 0:
            ntry = numgraphs * 10
            cands = _new_edges(graph, fix_outgroup)
            plus = [_insert_admix(graph, names, rng, rng.choice(cands)) for _ in range(ntry)] if cands else []
            plus = [g for g in plus if g is not None and allowed(g)]
            if not plus:
                continue
            scored = sorted(((fit(g).score, g) for g in plus), key=lambda x: x[0])
            score, g = scored[0]
            mut = "plusnadmix"
            if _numadmix(g) > max_admix:
                mut = "plusminusn"
                adm = [(p, a) for a in _admix_nodes(g) for p in _parents(g, a)]
                minus = [_delete_admix(g, names, *rng.choice(adm)) for _ in range(min(100, ntry))]
                minus = [m for m in minus if m is not None and allowed(m)]
                if not minus:
                    continue
                score, g = min(((fit(m).score, m) for m in minus), key=lambda x: x[0])
            if score > min(m.score for m in models.values()) * 0.95 or key(g) in models:
                gimp += 1
                continue
            new = [(g, mut)]
        else:
            neg = [e for e in fit(graph, constrained=False).edges if e.weight < 0]
            new = []
            if neg:
                k = numgraphs // 2
                cutoff = sorted(e.weight for e in neg)[min(k, len(neg)) - 1] if k else -float("inf")
                worst = [e for e in neg if e.weight <= cutoff]
                wfun = rng.choice([rearrange_negadmix3, replace_admix_with_random])
                for e in worst:
                    f = wfun if e.kind == "admix" else rearrange_negdrift
                    try:
                        g = f(graph, names, rng, e.from_node, e.to_node)
                    except (ValueError, IndexError):
                        g = None
                    if g is not None:
                        new.append((g, f.__name__))
            remaining = numgraphs - len(new)
            if remaining > 0 and neg:
                swaps = swap_negdrift(graph, names, rng, [(e.from_node, e.to_node) for e in neg])
                for g in rng.sample(swaps, min(len(swaps), -(-remaining // 2))):
                    new.append((g, "swap_negdrift"))
                remaining = numgraphs - len(new)
            for name in rng.choices(list(MUTATIONS), k=max(remaining, 0)):
                try:
                    g = MUTATIONS[name](graph, names, rng, fix_outgroup=fix_outgroup)
                except (ValueError, IndexError):
                    g = None
                if g is not None:
                    new.append((g, name))
            seen = set()
            kept = []
            for g, mut in new:
                if not _ok(g, pops) or not allowed(g):
                    continue
                if outgroup is not None and _outpop(g) != outgroup:
                    continue
                k = key(g)
                if k in models or k in seen:
                    continue
                seen.add(k)
                kept.append((g, mut))
            new = kept
            if not new:
                # as admixtools2: a generation without new graphs does not
                # count as one without improvement
                continue
        fitted = sorted((add(g, gen, mut) for g, mut in new), key=lambda m: m.score)
        if fitted[0].score < best:
            best = fitted[0].score
            gimp = 0
        else:
            gimp += 1
    return sorted(models.values(), key=lambda m: m.score), len(models)


def find_graphs(f2_blocks: _core.PyF2Blocks, populations: Sequence[str], nadmix: int = 0,
                outgroup: str | None = None, method: str = "stepwise", **kwargs: Any):
    """Search for the admixture graphs that best fit the data.

    Parameters
    ----------
    f2_blocks
        f2 blocks with every population.
    populations
        The graph's leaves; ``populations[0]`` is the f3 base of every fit.
    nadmix
        Number of admixture nodes.
    outgroup
        Only graphs with this population directly below the root.
    method
        ``"stepwise"`` (default) and ``"exhaustive"``: see
        ``popmix._core.find_graphs``. ``"evolutionary"``: admixtools2's
        ``find_graphs`` search (module docs), with its options
        ``initgraph``, ``numgraphs``, ``stop_gen``, ``stop_gen2``,
        ``stop_score``, ``stop_sec``, ``plusminus_generations``,
        ``max_admix``, ``admix_constraints``, ``event_constraints``,
        ``reject_f4z`` and ``seed``.
    **kwargs
        Options of the method, and of qpgraph (``numstart``, ``diag``,
        ``diag_f3``).

    Returns
    -------
    tuple
        (graphs, n_evaluated): the graphs, best first, and the number of
        graphs fitted. The evolutionary search returns every graph it fitted.
    """
    if method == "evolutionary":
        return evolve_graphs(f2_blocks, populations, nadmix, outgroup, **kwargs)
    return _core.find_graphs(f2_blocks, list(populations), nadmix, outgroup, method, **kwargs)
