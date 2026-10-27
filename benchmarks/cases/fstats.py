"""Benchmark cases: convertf, f2 blocks, qpAdm (one model, rotation), qpdstat, resampling."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ._common import (
    CONVERTF,
    POPS,
    QPADM,
    Case,
    _convertf_par,
    _eigenstrat_rows,
    _f2_mismatches,
    _py,
    _r,
)
from .graphs import GRAPH_POPS, _graph_file


def _agree_convertf(ours: Path, theirs: Path) -> dict:
    a, b = _eigenstrat_rows(ours / "out"), _eigenstrat_rows(theirs / "out")
    differ = sum(a.get(k) != v for k, v in b.items()) + len(a.keys() - b.keys())
    return {"snps": len(b), "differing_snps": differ, "pass": differ == 0}


CONVERTF_CASE = Case(
    name="convertf", analysis="packed to EIGENSTRAT text", reference="AdmixTools 8.0.2 convertf",
    popmix=lambda i, o, t: ["popmix", "--threads", str(t), "convertf", "-i", str(i), "-o", str(o / "out"),
                            "-of", "eigenstrat"],
    ref=lambda i, o: [CONVERTF, "-p", str(_convertf_par(i, o / "out", "EIGENSTRAT"))],
    agree=_agree_convertf,
)


_R_F2 = """
suppressMessages(library(admixtools))
a <- commandArgs(TRUE)
pops <- strsplit(a[2], ",")[[1]]
f2 <- f2_from_geno(a[1], pops = pops, verbose = FALSE)
writeBin(as.vector(f2), file.path(a[3], "f2.bin"))
jsonlite::write_json(list(dim = dim(f2), pops = dimnames(f2)[[1]]), file.path(a[3], "f2.json"),
                     auto_unbox = TRUE)
"""


def _agree_f2(ours: Path, theirs: Path) -> dict:
    import popmix as pm

    meta = json.loads((theirs / "f2.json").read_text())
    n, _, nb = meta["dim"]
    ref = np.fromfile(theirs / "f2.bin", dtype="<f8").reshape((nb, n, n)).transpose(1, 2, 0)
    f2 = pm.read_f2(str(ours / "f2"))
    worst_abs, differing = 0.0, 0
    for i, a in enumerate(meta["pops"]):
        for j, b in enumerate(meta["pops"]):
            if i < j:
                got, want = np.asarray(f2.f2_per_block(a, b)), ref[i, j]
                differing += _f2_mismatches(got, want)
                worst_abs = max(worst_abs, float(np.nanmax(np.abs(got - want), initial=0.0)))
    ok = list(f2.block_lengths) and len(f2.block_lengths) == nb
    return {"blocks": nb, "max_abs": worst_abs, "differing": differing,
            "pass": bool(ok) and differing == 0}


F2_CASE = Case(
    name="extract_f2", analysis=f"f2 blocks, {len(POPS)} populations", reference="admixtools2 2.0.10 f2_from_geno",
    popmix=lambda i, o, t: ["popmix", "--threads", str(t), "extract-f2", "-i", str(i), "-o", str(o / "f2"),
                            "--populations", ",".join(POPS), "--pseudohaploid", "admixtools2"],
    ref=lambda i, o: _r(_R_F2, o, str(i), ",".join(POPS)),
    agree=_agree_f2,
)


_R_QPADM = """
suppressMessages(library(admixtools))
a <- commandArgs(TRUE)
left <- strsplit(a[3], ",")[[1]]; right <- strsplit(a[4], ",")[[1]]
# f2 blocks first, as popmix qpadm --input does; qpadm(prefix) would compute
# each f4 from the genotypes on its own SNPs (DEVIATIONS.md)
f2 <- f2_from_geno(a[1], pops = c(a[2], left, right), verbose = FALSE)
r <- qpadm(f2, left = left, right = right, target = a[2], verbose = FALSE)
jsonlite::write_json(list(weights = r$weights$weight, se = r$weights$se,
                          p = r$rankdrop$p[1], chisq = r$rankdrop$chisq[1]),
                     file.path(a[5], "qpadm.json"), digits = NA)
"""


def _agree_qpadm(ours: Path, theirs: Path) -> dict:
    got = json.loads((ours / "qpadm.json").read_text())
    want = json.loads((theirs / "qpadm.json").read_text())
    w = np.array([got["weights"][s] for s in got["sources"]])
    se = np.array([got["std_errors"][s] for s in got["sources"]])
    rel = lambda x, y: float(np.max(np.abs(np.asarray(x) - np.asarray(y)) / np.abs(np.asarray(y))))  # noqa: E731
    out = {"weights_rel": rel(w, want["weights"]), "se_rel": rel(se, want["se"]),
           "p_rel": rel([got["p_value"]], [want["p"]])}
    out["pass"] = out["weights_rel"] <= 1e-6 and out["se_rel"] <= 1e-6 and out["p_rel"] <= 1e-4
    return out


def _popmix_qpadm(i: Path, o: Path, t: int, *extra: str) -> list[str]:
    q = QPADM
    return ["sh", "-c", " ".join(["popmix", "--threads", str(t), "qpadm", *extra, "-i", str(i), "--target", q["target"],
                                  "--left", ",".join(q["left"]), "--right", ",".join(q["right"]),
                                  "--format", "json", ">", str(o / "qpadm.json")])]


QPADM_CASE = Case(
    name="qpadm", analysis="qpAdm from genotypes (2 sources, 7 right)", reference="admixtools2 2.0.10 qpadm",
    popmix=_popmix_qpadm,
    ref=lambda i, o: _r(_R_QPADM, o, str(i), QPADM["target"], ",".join(QPADM["left"]), ",".join(QPADM["right"])),
    agree=_agree_qpadm,
)


DSTAT = {"pop1": ["Chimp"], "pop2": ["Mbuti"], "pop3": ["French", "Han", "Papuan", "Karitiana", "Sardinian"],
         "pop4": ["Yoruba", "Basque", "Turkey_N", "Russia_Samara_EBA_Yamnaya"]}


_R_QPDSTAT = """
suppressMessages(library(admixtools))
a <- commandArgs(TRUE)
p <- lapply(a[2:5], function(x) strsplit(x, ",")[[1]])
r <- qpdstat(a[1], p[[1]], p[[2]], p[[3]], p[[4]], verbose = FALSE)
jsonlite::write_json(r, file.path(a[6], "qpdstat.json"), digits = NA)
"""


def _agree_qpdstat(ours: Path, theirs: Path) -> dict:
    key = ("pop1", "pop2", "pop3", "pop4")
    want = {tuple(r[k] for k in key): r for r in json.loads((theirs / "qpdstat.json").read_text())}
    got = {tuple(r[k] for k in key): r for r in json.loads((ours / "qpdstat.json").read_text())}
    est = max(abs(got[k]["est"] - w["est"]) / abs(w["est"]) for k, w in want.items())
    se = max(abs(got[k]["se"] - w["se"]) / abs(w["se"]) for k, w in want.items())
    return {"statistics": len(want), "est_rel": est, "se_rel": se,
            "pass": got.keys() == want.keys() and est <= 1e-10 and se <= 1e-6}


QPDSTAT_CASE = Case(
    name="qpdstat", analysis="20 f4-statistics from genotypes", reference="admixtools2 2.0.10 qpdstat",
    popmix=lambda i, o, t: ["sh", "-c", f"popmix --threads {t} qpdstat -i {i} "
                                        + " ".join(f"--{k} {','.join(v)}" for k, v in DSTAT.items())
                                        + f" --format json > {o / 'qpdstat.json'}"],
    ref=lambda i, o: _r(_R_QPDSTAT, o, str(i), *(",".join(v) for v in DSTAT.values())),
    agree=_agree_qpdstat,
)


ROTATE = {"target": "Spain_EN",
          "leftright": ["Turkey_N", "Serbia_IronGates_Mesolithic", "Russia_Samara_EBA_Yamnaya", "Iran_GanjDareh_N",
                        "Israel_Natufian"],
          "rightfix": ["Mbuti", "Papuan", "Karitiana", "Han"]}


_R_ROTATE = """
suppressMessages(library(admixtools))
future::plan("sequential")
a <- commandArgs(TRUE)
lr <- strsplit(a[3], ",")[[1]]; rf <- strsplit(a[4], ",")[[1]]
f2 <- f2_from_geno(a[1], pops = c(a[2], lr, rf), verbose = FALSE)
r <- qpadm_rotate(f2, lr, a[2], rf, verbose = FALSE)
jsonlite::write_json(r, file.path(a[5], "rotate.json"), digits = NA)
"""


_PY_ROTATE = """
import json, sys
import popmix as pm
t, prefix, target, lr, rf, out = sys.argv[1:]
pm.init_threads(int(t))
lr, rf = lr.split(","), rf.split(",")
f2 = pm.extract_f2(prefix, [target, *lr, *rf])
r = pm.qpadm_rotate(f2, lr, target, rf)
rows = [{k: (list(v) if k in ("left", "right") else v.item() if hasattr(v, "item") else v) for k, v in row.items()}
        for row in r.to_dict("records")]
open(out + "/rotate.json", "w").write(json.dumps(rows))
"""


def _agree_rotate(ours: Path, theirs: Path) -> dict:
    def key(r: dict) -> tuple:
        return tuple(r["left"]), tuple(r["right"])

    got = {key(r): r for r in json.loads((ours / "rotate.json").read_text())}
    want = {key(r): r for r in json.loads((theirs / "rotate.json").read_text())}
    if got.keys() != want.keys():
        return {"models": len(want), "pass": False}
    # tests/golden/test_qpadm_models_admixtools2.py
    exact = all((got[k]["f4rank"], got[k]["dof"], got[k]["feasible"]) == (w["f4rank"], w["dof"], w["feasible"])
                for k, w in want.items())
    chisq = max(abs(got[k]["chisq"] - w["chisq"]) / abs(w["chisq"]) for k, w in want.items())
    p = max(0.0 if got[k]["p"] == w["p"] else abs(got[k]["p"] - w["p"]) / abs(w["p"]) for k, w in want.items())
    return {"models": len(want), "rank_dof_feasible_equal": exact, "chisq_rel": chisq, "p_rel": p,
            "pass": exact and chisq <= 1e-6 and p <= 1e-4}


def _rotate_args() -> list[str]:
    return [ROTATE["target"], ",".join(ROTATE["leftright"]), ",".join(ROTATE["rightfix"])]


ROTATE_CASE = Case(
    name="qpadm_rotate", analysis=f"qpAdm rotation, {len(ROTATE['leftright'])} populations between left and right, "
                                  f"{len(ROTATE['rightfix'])} fixed right",
    reference="admixtools2 2.0.10 qpadm_rotate",
    popmix=lambda i, o, t: _py(_PY_ROTATE, o, t, str(i), *_rotate_args()),
    ref=lambda i, o: _r(_R_ROTATE, o, str(i), *_rotate_args()),
    agree=_agree_rotate,
)


RESAMPLE_BOOT, RESAMPLE_SEED = 10, 3


_R_RESAMPLE = """
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
edges <- as.matrix(read.table(a[2], col.names = c("from", "to")))
f2 <- f2_from_geno(a[1], pops = strsplit(a[3], ",")[[1]], verbose = FALSE)
# as tests/golden/admixtools2/generate_resample.R: the draws from set.seed(),
# qpgraph's random starts from seed = 1
set.seed(as.integer(a[5]))
r <- qpgraph_resample_snps(f2, boot = as.integer(a[4]), graph = edges, seed = 1)
jsonlite::write_json(list(score = r$score, sel = r$sel), file.path(a[6], "resample.json"), digits = NA)
"""


_PY_RESAMPLE = """
import json, sys
import popmix as pm
t, prefix, graph, pops, boot, seed, out = sys.argv[1:]
pm.init_threads(int(t))
edges = [tuple(line.split()) for line in open(graph) if line.strip()]
f2 = pm.extract_f2(prefix, pops.split(","))
r = pm.qpgraph_resample_snps(f2, edges, boot=int(boot), boot_seed=int(seed), seed=1)
open(out + "/resample.json", "w").write(json.dumps({"score": [float(s) for s in r["score"]],
                                                    "sel": [[int(b) + 1 for b in s] for s in r["sel"]]}))
"""


def _agree_resample(ours: Path, theirs: Path) -> dict:
    got = json.loads((ours / "resample.json").read_text())
    want = json.loads((theirs / "resample.json").read_text())
    same_draws = got["sel"] == want["sel"]
    rel = max(abs(a - b) / abs(b) for a, b in zip(got["score"], want["score"]))  # noqa: B905
    # tests/golden/test_resample_admixtools2.py: optimized scores rtol 1e-8
    return {"replicates": len(want["score"]), "same_blocks_drawn": same_draws, "score_rel": rel,
            "pass": same_draws and len(got["score"]) == len(want["score"]) and rel <= 1e-8}


def _resample_args(o: Path) -> list[str]:
    return [str(_graph_file(o)), ",".join(GRAPH_POPS), str(RESAMPLE_BOOT), str(RESAMPLE_SEED)]


RESAMPLE_CASE = Case(
    name="qpgraph_resample", analysis=f"qpGraph bootstrap, {RESAMPLE_BOOT} resamples of the blocks, "
                                      f"{len(GRAPH_POPS)} populations, 1 admixture",
    reference="admixtools2 2.0.10 qpgraph_resample_snps",
    popmix=lambda i, o, t: _py(_PY_RESAMPLE, o, t, str(i), *_resample_args(o)),
    ref=lambda i, o: _r(_R_RESAMPLE, o, str(i), *_resample_args(o)),
    agree=_agree_resample,
)
