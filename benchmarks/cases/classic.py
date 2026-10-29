"""Benchmark cases: the classic AdmixTools programs (qpDstat, qp3Pop, qpAdm, qpGraph)."""

from __future__ import annotations

import json
import re
from pathlib import Path

from ._common import CLASSIC, QPADM, Case
from .fstats import DSTAT, _popmix_qpadm
from .graphs import GRAPH, GRAPH_POPS, _graph_file

CLASSIC_BIN = CLASSIC / "bin"


PRINTED6, PRINTED3 = 5e-7 + 1e-12, 5e-4 + 1e-9


def _classic_par(i: Path, o: Path, name: str, **keys) -> Path:
    par = o / f"{name}.par"
    par.write_text("".join(f"{k}: {v}\n" for k, v in {"genotypename": f"{i}.geno", "snpname": f"{i}.snp",
                                                       "indivname": f"{i}.ind", **keys}.items()), encoding="utf-8")
    return par


def _classic(prog: str, par: Path, *args: str) -> list[str]:
    # the classic programs print their results to stdout
    return ["sh", "-c", f"{CLASSIC_BIN / prog} -p {par} {' '.join(args)} > {par.with_suffix('.out')}"]


def _results(out: Path) -> list[list[str]]:
    return [line.split()[1:] for line in out.read_text(encoding="utf-8").splitlines() if line.split()[:1] == ["result:"]]


def _quads() -> list[tuple[str, ...]]:
    return [(a, b, c, d) for a in DSTAT["pop1"] for b in DSTAT["pop2"] for c in DSTAT["pop3"] for d in DSTAT["pop4"]]


def _classic_qpdstat(i: Path, o: Path) -> list[str]:
    (o / "quads").write_text("".join(" ".join(q) + "\n" for q in _quads()), encoding="utf-8")
    return _classic("qpDstat", _classic_par(i, o, "qpdstat", popfilename=o / "quads", f4mode="YES"))


def _agree_classic_qpdstat(ours: Path, theirs: Path) -> dict:
    got = {tuple(r[k] for k in ("pop1", "pop2", "pop3", "pop4")): r
           for r in json.loads((ours / "qpdstat.json").read_text(encoding="utf-8"))}
    # result: pop1 pop2 pop3 pop4 f4 Z BABA ABBA n
    want = {tuple(r[:4]): (float(r[4]), float(r[5]), int(r[8])) for r in _results(theirs / "qpdstat.out")}
    if got.keys() != want.keys():
        return {"statistics": len(want), "pass": False}
    est = max(abs(got[k]["est"] - w[0]) for k, w in want.items())
    z = max(abs(got[k]["z"] - w[1]) for k, w in want.items())
    n = all(got[k]["n"] == w[2] for k, w in want.items())
    return {"statistics": len(want), "est_abs": est, "z_abs": z, "snps_equal": n,
            "pass": est <= PRINTED6 and z <= PRINTED3 and n}


CLASSIC_QPDSTAT_CASE = Case(
    name="classic_qpdstat", analysis="20 f4-statistics from genotypes", reference="ADMIXTOOLS 7.0.2 qpDstat",
    popmix=lambda i, o, t: ["sh", "-c", f"popmix --threads {t} qpdstat -i {i} --pseudohaploid diploid "
                                        + " ".join(f"--{k} {','.join(v)}" for k, v in DSTAT.items())
                                        + f" --format json > {o / 'qpdstat.json'}"],
    ref=_classic_qpdstat, agree=_agree_classic_qpdstat,
)


# outgroup f3(Mbuti; A, B)
QP3POP = {"pop1": ["Mbuti"], "pop2": ["French", "Han", "Papuan", "Karitiana"],
          "pop3": ["Basque", "Sardinian", "Turkey_N", "Russia_Samara_EBA_Yamnaya"]}


def _classic_qp3pop(i: Path, o: Path) -> list[str]:
    # qp3Pop's lines are "source1 source2 target"
    (o / "triples").write_text("".join(f"{b} {c} {a}\n" for a in QP3POP["pop1"] for b in QP3POP["pop2"]
                                       for c in QP3POP["pop3"]), encoding="utf-8")
    return _classic("qp3Pop", _classic_par(i, o, "qp3pop", popfilename=o / "triples", outgroupmode="YES"))


def _agree_classic_qp3pop(ours: Path, theirs: Path) -> dict:
    got = {(r["pop1"], r["pop2"], r["pop3"]): r for r in json.loads((ours / "qp3pop.json").read_text(encoding="utf-8"))}
    # result: source1 source2 target f3 se Z n; outgroupmode prints f3 and se x 1000
    want = {(r[2], r[0], r[1]): (float(r[3]) / 1000, float(r[4]) / 1000, float(r[5])) for r in
            _results(theirs / "qp3pop.out")}
    if got.keys() != want.keys():
        return {"statistics": len(want), "pass": False}
    est = max(abs(got[k]["est"] - w[0]) for k, w in want.items())
    se = max(abs(got[k]["se"] - w[1]) for k, w in want.items())
    z = max(abs(got[k]["z"] - w[2]) for k, w in want.items())
    # qp3Pop counts its SNPs otherwise, so n is not compared
    return {"statistics": len(want), "est_abs": est, "se_abs": se, "z_abs": z,
            "pass": est <= PRINTED6 / 1000 and se <= PRINTED6 / 1000 and z <= PRINTED3}


CLASSIC_QP3POP_CASE = Case(
    name="classic_qp3pop", analysis="16 outgroup f3-statistics from genotypes", reference="ADMIXTOOLS 7.0.2 qp3Pop",
    popmix=lambda i, o, t: ["sh", "-c", f"popmix --threads {t} qp3pop -i {i} --pseudohaploid diploid --outgroup-mode "
                                        + " ".join(f"--{k} {','.join(v)}" for k, v in QP3POP.items())
                                        + f" --format json > {o / 'qp3pop.json'}"],
    ref=_classic_qp3pop, agree=_agree_classic_qp3pop,
)


def _classic_qpadm(i: Path, o: Path) -> list[str]:
    q = QPADM
    (o / "left").write_text("".join(f"{p}\n" for p in [q["target"], *q["left"]]), encoding="utf-8")
    (o / "right").write_text("".join(f"{p}\n" for p in q["right"]), encoding="utf-8")
    return _classic("qpAdm", _classic_par(i, o, "qpadm", popleft=o / "left", popright=o / "right"))


def _agree_classic_qpadm(ours: Path, theirs: Path) -> dict:
    got = json.loads((ours / "qpadm.json").read_text(encoding="utf-8"))
    text = (theirs / "qpadm.out").read_text(encoding="utf-8")
    weights = [float(x) for x in re.search(r"^best coefficients:(.*)$", text, re.M).group(1).split()]
    se = [float(x) for x in re.search(r"^\s*std\. errors:(.*)$", text, re.M).group(1).split()]
    rank = len(QPADM["left"]) - 1
    m = re.search(rf"^f4rank:\s+{rank} dof:\s+(\d+) chisq:\s+(\S+) tail:\s+(\S+)", text, re.M)
    chisq, p = float(m.group(2)), float(m.group(3))
    w = [got["weights"][s] for s in got["sources"]]
    s = [got["std_errors"][x] for x in got["sources"]]
    out = {"weights_abs": max(abs(a - b) for a, b in zip(w, weights)),  # noqa: B905
           "se_abs": max(abs(a - b) for a, b in zip(s, se)),  # noqa: B905
           "chisq_abs": abs(got["chi_sq"] - chisq), "p_rel": abs(got["p_value"] - p) / p}
    # qpAdm keeps another SNP set from admixtools2's and popmix's (803,891
    # SNPs here, against 675,980), so on this model the two differ by more
    # than the printed precision (weights by 1e-3, chi-squared by 1.5%): the
    # case records how far. Classic qpAdm keeps a SNP when every basis
    # f4(L0, Li; R0, Rj) can be evaluated, an f4 with a factor exactly 0
    # counting as 0 even when the other factor's population has no data, and
    # keeps monomorphic SNPs; on the SNPs every population has, its standard
    # errors and chi-squared still differ from admixtools2's (DEVIATIONS.md,
    # "ADMIXTOOLS 7.0.2 (qpAdm)"). No popmix option reproduces it.
    # weights, SEs and chi-squared printed to 3 decimals; the tail probability
    # to 10 significant digits, so CLAUDE.md's p-value tolerance applies
    out["pass"] = (out["weights_abs"] <= PRINTED3 and out["se_abs"] <= PRINTED3 and out["chisq_abs"] <= PRINTED3
                   and out["p_rel"] <= 1e-4)
    return out


CLASSIC_QPADM_CASE = Case(
    name="classic_qpadm", analysis="qpAdm from genotypes (2 sources, 7 right)", reference="ADMIXTOOLS 7.0.2 qpAdm",
    popmix=lambda i, o, t: _popmix_qpadm(i, o, t, "--pseudohaploid", "diploid"),
    ref=_classic_qpadm, agree=_agree_classic_qpadm,
)


def _classic_graph(o: Path) -> Path:
    """GRAPH in qpGraph's format: drift edges, and an admix line for each
    node with two parents."""
    parents: dict[str, list[str]] = {}
    for a, b in GRAPH:
        parents.setdefault(b, []).append(a)
    lines = [f"root {GRAPH[0][0]}", *(f"label {p} {p}" for p in GRAPH_POPS)]
    for k, (a, b) in enumerate(GRAPH):
        if len(parents[b]) == 1:
            lines.append(f"edge e{k} {a} {b}")
    lines += [f"admix {b} {' '.join(p)}" for b, p in parents.items() if len(p) == 2]
    path = o / "graph.classic"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _classic_qpgraph(i: Path, o: Path) -> list[str]:
    # tests/golden/admixtools_classic (removed in 9ed980b) ran qpGraph so
    par = _classic_par(i, o, "qpgraph", outpop="NULL", useallsnps="NO", blgsize=0.05, forcezmode="YES",
                       lsqmode="NO", diag=".0001", bigiter=6, hires="YES", lambdascale=1)
    return _classic("qpGraph", par, "-g", str(_classic_graph(o)), "-o", str(o / "qpgraph.ggg"))


def _agree_classic_qpgraph(ours: Path, theirs: Path) -> dict:
    got = json.loads((ours / "qpgraph.json").read_text(encoding="utf-8"))["score"]
    want = float(re.search(r"final score:\s+(\S+)", (theirs / "qpgraph.out").read_text(encoding="utf-8")).group(1))
    # printed to 3 decimals; both optimize, R's L-BFGS-B to about 1e-8 of the score
    return {"score": want, "score_abs": abs(got - want), "pass": abs(got - want) <= PRINTED3 + 1e-6 * abs(want)}


CLASSIC_QPGRAPH_CASE = Case(
    name="classic_qpgraph", analysis="qpGraph, 6 populations, 1 admixture", reference="ADMIXTOOLS 7.0.2 qpGraph",
    popmix=lambda i, o, t: ["sh", "-c", f"popmix --threads {t} qpgraph -i {i} --graph {_graph_file(o)} "
                                        f"--pseudohaploid diploid --numstart 10 --seed 1 --format json "
                                        f"> {o / 'qpgraph.json'}"],
    ref=_classic_qpgraph, agree=_agree_classic_qpgraph,
)
