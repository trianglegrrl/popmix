"""Admixture graphs in Graphviz DOT format, as admixtools2's ``write_dot``
and ``parse_dot``.

Drift edges are labelled with their length times 1000 and admixture edges,
dotted, with their weight as a percentage; internal nodes are points.

admixtools2 removes "." and "-" from the node names on edges (Han.DG
becomes HanDG) but not from the node declarations, so a dotted internal
name is declared under one name and used under another (DEVIATIONS.md,
ADMIXTOOLS2-0004). popmix quotes names and keeps them; ``compat=
"admixtools2-2.0.10"`` writes admixtools2's text exactly.
"""

from __future__ import annotations

import math
import os
import re
from typing import Any

COMPAT = "admixtools2-2.0.10"


def _edges(graph: Any) -> list[tuple[str, str, str, float]]:
    """(from, to, type, weight) for a qpgraph result or an edge list."""
    if hasattr(graph, "edges") and not isinstance(graph, (list, tuple)):
        return [(e.from_node, e.to_node, e.kind, e.weight) for e in graph.edges]
    rows = [tuple(e) for e in graph]
    if all(len(r) >= 4 for r in rows):
        return [(a, b, t, w) for a, b, t, w, *_ in rows]
    indeg: dict[str, int] = {}
    for _, b, *_ in rows:
        indeg[b] = indeg.get(b, 0) + 1
    return [(a, b, "admix" if indeg[b] > 1 else "edge", 0.0) for a, b, *_ in rows]


def _num(x: float) -> str:
    """R's as.character for a double."""
    return f"{x:.15g}"


def _rounded(x: float) -> str:
    return "NA" if x is None or math.isnan(x) else str(int(round(x)))


def write_dot(graph: Any, path: str | os.PathLike | None = None, *, title: str = "",
              size: tuple[float, float] = (7.5, 10), compat: str | None = None) -> str:
    """The graph in DOT format (admixtools2 ``write_dot``).

    Parameters
    ----------
    graph
        A ``qpgraph`` result (edges labelled with their fit), or the graph
        as (from, to) pairs (unlabelled, as admixtools2 writes an igraph), or
        as (from, to, type, weight) rows.
    path
        Also write it to this file.
    title, size
        The graph's label and Graphviz page size in inches.
    compat
        ``"admixtools2-2.0.10"``: admixtools2's text exactly.

    Returns
    -------
    str
        The DOT text.
    """
    if compat not in (None, COMPAT):
        raise ValueError(f"unknown compat {compat!r}; expected None or {COMPAT!r}")
    edges = _edges(graph)
    froms = [e[0] for e in edges]
    tos = [e[1] for e in edges]
    leaves = set(tos) - set(froms)
    internal = [v for v in dict.fromkeys(froms + tos) if v not in leaves]

    if compat:
        def name(v: str) -> str:
            return re.sub(r"[.-]", "", v)

        def decl(v: str) -> str:
            return v
        arrow, sep = "  ->  ", " "
    else:
        def name(v: str) -> str:
            return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
        decl = name
        arrow, sep = " -> ", " "

    lines = [f'digraph G {{\nlabel = "{title}";\nlabelloc=t;\nlabeljust=l;',
             f'size = "{_num(size[0])},{_num(size[1])}";',
             "\n".join(f"{decl(v)} [shape = point];" for v in internal)]
    body = []
    for a, b, kind, w in edges:
        if kind == "admix":
            lab = f' [ style=dotted, label = "{_rounded(w * 100)}%" ];'
        else:
            lab = f' [ label = "{_rounded(w * 1000)}" ];'
        body.append(f"{name(a)}{arrow}{name(b)}{sep}{lab}")
    text = "\n".join(lines) + "\n" + "\n".join(body) + "\n}\n"
    if path is not None:
        with open(path, "w") as fh:
            fh.write(text)
    return text


def parse_dot(source: str | os.PathLike) -> list[tuple[str, str]]:
    """The edges of a DOT file or text (admixtools2 ``parse_dot``): each
    "a -> b" line, names unquoted.

    Parameters
    ----------
    source
        Path to a DOT file, or DOT text (a string containing "->" is read as
        text).

    Returns
    -------
    list of tuple of str
        (from, to) for each edge line, in file order; labels and other
        attributes are dropped.
    """
    text = source
    if isinstance(source, os.PathLike) or (isinstance(source, str) and "->" not in source and os.path.exists(source)):
        with open(source) as fh:
            text = fh.read()
    out = []
    for line in str(text).splitlines():
        if "->" not in line:
            continue
        a, b = line.split("->", 1)
        out.append((_token(a), _token(b)))
    return out


def _token(s: str) -> str:
    s = s.strip()
    if s.startswith('"'):
        m = re.match(r'"((?:[^"\\]|\\.)*)"', s)
        if m:
            return re.sub(r"\\(.)", r"\1", m.group(1))
    return s.split()[0]
