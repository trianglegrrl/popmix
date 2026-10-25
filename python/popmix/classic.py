"""Read the output of classic ADMIXTOOLS programs, as admixtools2's
``parse_qpdstat_output``, ``parse_qp3pop_output``, ``parse_qpadm_output``,
``parse_qpf4ratio_output``, ``parse_qpgraph_output``,
``parse_qpgraph_output_edges``, ``parse_qpgraph_graphfile`` and
``parse_qpfstats_output`` do: the same tables, with the same derived
columns (``se`` from the statistic and its Z-score, two-sided normal
``p.value``).

``parse_qpgraph_graphfile`` returns the populations' names where admixtools2
returns the graph's node labels (DEVIATIONS.md, ADMIXTOOLS2-0005);
``compat="admixtools2-2.0.10"`` returns the labels.
"""

from __future__ import annotations

import math
import os
import re
from typing import Any

import pandas as pd

COMPAT = "admixtools2-2.0.10"


def _lines(path: str | os.PathLike) -> list[str]:
    with open(path) as fh:
        return fh.read().splitlines()


def _squish(s: str) -> str:
    return " ".join(s.split())


def _num(x: str) -> Any:
    try:
        v = float(x)
    except ValueError:
        return x
    return int(v) if re.fullmatch(r"-?\d+", x) else v


def _ztop(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2))


def parse_qpdstat_output(path: str | os.PathLike) -> pd.DataFrame:
    """qpDstat's ``result:`` lines: pop1-pop4, f4 (or D), Z, [best,] BABA,
    ABBA, numsnps, se and p.value."""
    rows = [_squish(x).split(" ") for x in _lines(path) if x.startswith("result: ")]
    if not rows:
        raise ValueError(f"{path}: no result lines")
    nam1 = ["result", "pop1", "pop2", "pop3", "pop4", "f4", "Z"]
    nam2 = ["BABA", "ABBA", "numsnps"]
    names = {10: nam1 + nam2, 11: nam1 + ["best"] + nam2}.get(len(rows[0]))
    if names is None:
        raise ValueError(f"{path}: unexpected number of columns in output")
    t = pd.DataFrame([[_num(v) for v in r] for r in rows], columns=names).drop(columns="result")
    t["se"] = t.f4 / t.Z
    t["p.value"] = [_ztop(z) for z in t.Z]
    return t


def parse_qp3pop_output(path: str | os.PathLike) -> pd.DataFrame:
    """qp3Pop's ``result:`` lines: source1, source2, target, f3, se (as
    f3 / Z), Z, numsnps and p.value."""
    lines = [x.replace(" no data", "result:") for x in _lines(path)]
    rows = [_squish(x).split(" ") for x in lines if x.startswith(" result: ")]
    if not rows:
        raise ValueError(f"{path}: no result lines")
    names = ["result", "source1", "source2", "target", "f3", "se", "Z", "numsnps"]
    t = pd.DataFrame([[_num(v) for v in r] for r in rows], columns=names).drop(columns="result")
    t["se"] = t.f3 / t.Z
    t["p.value"] = [_ztop(z) for z in t.Z]
    return t


def _first(lines: list[str], pattern: str, start: int = 0) -> int:
    for i in range(start, len(lines)):
        if re.search(pattern, lines[i]):
            return i
    raise ValueError(f"no line matching {pattern!r}")


def parse_qpadm_output(path: str | os.PathLike) -> dict[str, pd.DataFrame]:
    """qpAdm's weights (target, left, weight, mean, se, z) and popdrop
    (pat, wt, dof, chisq, p, one column per source, feasible)."""
    dat = _lines(path)
    lstart = _first(dat, "left pops:") + 1
    lend = _first(dat, r"^$", lstart + 1) - 1
    coefstart = _first(dat, r"^best coefficients:")
    sigstart = _first(dat, "fixed pat") + 1
    sigend = _first(dat, r"^best pat:") - 1
    target = dat[lstart]
    left = dat[lstart:lend + 1][1:]
    coefs = {}
    for name, line in zip(["weight", "mean", "se"], dat[coefstart:coefstart + 3]):  # noqa: B905
        parts = re.split(r" +", line)
        coefs[name] = [float(v) for v in parts[-(len(left) + 1):-1]]
    weights = pd.DataFrame({"target": target, "left": left, **coefs})
    weights["z"] = weights["mean"] / weights["se"]
    rows = [re.split(r" +", x)[1:] for x in dat[sigstart:sigend + 1]]
    cols = ["pat", "wt", "dof", "chisq", "p"] + left + ["feasible"]
    popdrop = pd.DataFrame([r[:len(cols)] + [""] * (len(cols) - len(r)) for r in rows], columns=cols)
    popdrop["feasible"] = popdrop.feasible != "infeasible"
    for c in cols:
        if c not in ("pat", "feasible"):
            popdrop[c] = popdrop[c].astype(float)
    return {"weights": weights, "popdrop": popdrop}


def parse_qpf4ratio_output(path: str | os.PathLike) -> pd.DataFrame:
    """qpF4ratio's ``result:`` lines: pop1-pop5, alpha, se, z and n."""
    rows = []
    for line in _lines(path):
        if " result: " not in line:
            continue
        p = re.split(r" +", line)
        rows.append([p[i - 1] for i in (3, 4, 5, 6, 10, 12, 13, 14, 15)])
    t = pd.DataFrame(rows, columns=["pop1", "pop2", "pop3", "pop4", "pop5", "alpha", "se", "z", "n"])
    for c in ("alpha", "se", "z"):
        t[c] = t[c].astype(float)
    t["n"] = t["n"].astype(int)
    return t


def parse_qpgraph_output_edges(path: str | os.PathLike, uselabels: bool = False) -> pd.DataFrame:
    """The fitted edges of a qpGraph graph output (``-o``): from, to, type
    (edge or admix) and weight, with population names (labels with
    ``uselabels``)."""
    dat = [_squish(x) for x in _lines(path)]
    nammap = {}
    for x in dat:
        if x.startswith("label"):
            p = x.split(" ")
            nammap[p[1]] = p[2]
    edges, admix1, admix2 = [], [], []
    for x in dat:
        p = x.split(" ")
        if not re.match(r"^(edge|ledge|redge|admix)", x):
            continue
        p += [None] * (6 - len(p))
        typ, name, frm, to, w, w2 = p[:6]
        if typ == "admix":
            admix1.append((frm, name, "admix", float(w)))
            admix2.append((to, name, "admix", float(w2)))
        else:
            edges.append((frm, to, "edge", float(w)))
    t = pd.DataFrame(edges + admix1 + admix2, columns=["from", "to", "type", "weight"])
    if not uselabels:
        t["from"] = t["from"].map(lambda v: nammap.get(v, v))
        t["to"] = t["to"].map(lambda v: nammap.get(v, v))
    return t


def parse_qpgraph_output(outfile: str | os.PathLike, logfile: str | os.PathLike | None = None) -> dict[str, Any]:
    """qpGraph's fit: ``edges`` (from the graph output, ``-o``), ``score``,
    ``f2`` (fit, est, diff, se, z), ``f3`` (from the ff3 tables) and ``f4``
    (the outliers, or None when some are missing, as admixtools2) from its
    log (standard output)."""
    edges = parse_qpgraph_output_edges(outfile)
    dat = [_squish(x) for x in _lines(logfile if logfile is not None else outfile)]
    score = next(float(x.split(" ")[2]) for x in dat if x.startswith("final score"))
    pops = [x.split(" ")[2] for x in dat if x.startswith("population:")]
    n = len(pops)
    f2rows = [x.split(" ") for x in dat if " f2: " in f" {x} " and x.split(" ")[2:3] == ["f2:"]]
    pairs = [(pops[i], pops[j]) for i in range(n - 1) for j in range(i + 1, n)]
    f2 = pd.DataFrame([[a, b] + [float(v) for v in r[3:8]] for (a, b), r in zip(pairs, f2rows)],  # noqa: B905
                      columns=["pop1", "pop2", "fit", "est", "diff", "se", "z"])

    def block(start: int, end: int) -> list[float]:
        return [float(v) for x in dat[start:end + 1] for v in x[2:].split(" ")]

    est0 = _first(dat, r"^ff3:") + 2
    est1 = _first(dat, r"^ff3sig\*10:") - 2
    sig0 = est1 + 4
    fit0 = _first(dat, r"^ff3fit:") + 2
    m = est1 - est0
    est, sig, fit = block(est0, est1), block(sig0, sig0 + m), block(fit0, fit0 + m)
    grid = [(a, b) for a in pops for b in pops]
    f3 = pd.DataFrame({"pop1": pops[0], "pop2": [a for a, _ in grid], "pop3": [b for _, b in grid],
                       "est": [v / 1000 for v in est], "fit": [v / 1000 for v in fit],
                       "se": [v / 1000 / 10 for v in sig]})
    f3["diff"] = f3.est - f3.fit
    # R: 0/0 is NaN, x/0 is +-Inf
    f3["z"] = [d / s if s else (math.nan if d == 0 else math.copysign(math.inf, d))
               for d, s in zip(f3["diff"], f3.se)]  # noqa: B905
    o0 = _first(dat, r"^outliers:") + 2
    o1 = _first(dat, r"^worst f-stat:") - 3
    nall = math.comb(math.comb(n, 2) + 1, 2)
    f4 = None
    if o1 - o0 + 1 == nall:
        pp = [(a, b) for a, b in pairs]
        quads = [pp[i] + pp[j] for i in range(len(pp)) for j in range(i, len(pp))]
        rows = []
        for q, x in zip(quads, dat[o0:o1 + 1]):  # noqa: B905
            v = [float(s) for s in x.split(" ")[4:9]]
            rows.append(list(q) + [v[0], v[1], -v[2], v[3], -v[4]])
        f4 = pd.DataFrame(rows, columns=["pop1", "pop2", "pop3", "pop4", "fit", "est", "diff", "se", "z"])
    return {"edges": edges, "score": score, "f2": f2, "f3": f3, "f4": f4}


def parse_qpgraph_graphfile(path: str | os.PathLike, uselabels: bool = False,
                            compat: str | None = None) -> pd.DataFrame:
    """A graph in qpGraph's format: from, to, and lower and upper bounds of
    locked admixture weights when there are any.

    A graph file names its nodes by labels (``label <population> <label>``);
    without ``uselabels`` they are replaced by the populations' names.
    ``compat="admixtools2-2.0.10"`` keeps the labels, as admixtools2 does."""
    if compat not in (None, COMPAT):
        raise ValueError(f"unknown compat {compat!r}; expected None or {COMPAT!r}")
    lines = _lines(path)
    label_to_pop = {}
    rows, locks = [], set()
    for line in lines:
        p = line.split()
        if not p:
            continue
        if p[0] == "label" and len(p) >= 3:
            label_to_pop[p[2]] = p[1]
        elif p[0] == "lock" and len(p) >= 2:
            locks.add(p[1])
        elif p[0] in ("edge", "ledge", "redge", "admix"):
            rows.append(p + [None] * (6 - len(p)))
    out = []
    for _, _, frm, to, *_ in (r for r in rows if r[0] != "admix"):
        out.append((frm, to, math.nan, math.nan))
    adm = [r for r in rows if r[0] == "admix"]
    # an admixture line gives both incoming edges of its node, in this order
    for _, name, frm, _, l1, _ in adm:
        lo = float(l1) / 100 if name in locks else math.nan
        out.append((frm, name, lo, lo))
    for _, name, _, to, _, l2 in adm:
        lo = float(l2) / 100 if name in locks else math.nan
        out.append((to, name, lo, lo))
    t = pd.DataFrame(out, columns=["from", "to", "lower", "upper"])
    if not uselabels and not compat:
        t["from"] = t["from"].map(lambda v: label_to_pop.get(v, v))
        t["to"] = t["to"].map(lambda v: label_to_pop.get(v, v))
    if t.lower.isna().all() and t.upper.isna().all():
        t = t.drop(columns=["lower", "upper"])
    return t


def parse_qpfstats_output(path: str | os.PathLike) -> dict[str, pd.DataFrame | None]:
    """qpfstats' ``jest``, ``fstat`` and ``basis`` tables, populations by
    name (None where the output has no such lines)."""
    dat = _lines(path)
    pops = [_squish(x).split(" ")[1] for x in dat if x.startswith("pop: ")]

    def rows(prefix: str, drop: list[int]) -> list[list[str]]:
        out = []
        for x in dat:
            if x.startswith(prefix):
                p = _squish(x).split(" ")
                out.append([v for i, v in enumerate(p) if i not in drop])
        return out

    def table(r: list[list[str]], names: list[str], popcols: list[str]) -> pd.DataFrame:
        t = pd.DataFrame([[float(v) for v in row] for row in r], columns=names)
        for c in popcols:
            t[c] = [pops[int(v)] for v in t[c]]
        return t

    jest = rows("jest", [0, 1, 2])
    fstat = rows("fstat:", [0])
    basis = rows("basis:", [0, 5])
    out = {
        "jest": table(jest, ["p11", "p12", "p21", "p22", "f21", "f22", "se"], ["p11", "p12", "p21", "p22"])
        if jest else None,
        "fstat": table(fstat, ["pop1", "pop2", "pop3", "pop4", "est", "se", "z"], ["pop1", "pop2", "pop3", "pop4"])
        if fstat else None,
        "basis": table(basis, ["pop1", "pop2", "f21", "se1", "f22", "se2"], ["pop1", "pop2"]),
    }
    return out
