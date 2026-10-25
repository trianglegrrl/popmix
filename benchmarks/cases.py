"""Benchmark cases: popmix and a reference tool on identical inputs.

Every case runs the popmix command line (one process, as a user would) and
the reference tool's own program, each writing its result to its own
directory; `agree` then compares the two results at the CLAUDE.md
tolerances. Inputs are made once, by the reference tools, in the cache.
"""

from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

HOME = Path.home()
CACHE = Path(os.environ.get("POPMIX_REFERENCE_CACHE", HOME / ".cache/popmix-references")) / "benchmarks"
AADR = Path("/references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB")
CONVERTF = os.environ.get("POPMIX_CONVERTF", str(HOME / ".local/share/popmix-references/admixtools-8.0.2/bin/convertf"))
RSCRIPT = os.environ.get("POPMIX_RSCRIPT", "/usr/bin/Rscript")

# 16 populations, present-day and ancient, ~450 individuals of AADR v66.p1
POPS = ["Mbuti", "Yoruba", "French", "Han", "Papuan", "Karitiana", "Basque", "Sardinian",
        "Russia_Samara_EBA_Yamnaya", "Turkey_N", "Serbia_IronGates_Mesolithic", "Iran_GanjDareh_N",
        "Israel_Natufian", "Morocco_Iberomaurusian", "Spain_EN", "Chimp"]
QPADM = {"target": "Spain_EN", "left": ["Turkey_N", "Serbia_IronGates_Mesolithic"],
         "right": ["Mbuti", "Papuan", "Karitiana", "Han", "Iran_GanjDareh_N", "Israel_Natufian",
                   "Morocco_Iberomaurusian"]}


def _convertf_par(src: Path, out: Path, fmt: str, **extra: str) -> Path:
    par = out.parent / f"{out.name}.par"
    keys = {"genotypename": f"{src}.geno", "snpname": f"{src}.snp", "indivname": f"{src}.ind",
            "outputformat": fmt, "genooutfilename": f"{out}.geno", "snpoutfilename": f"{out}.snp",
            "indoutfilename": f"{out}.ind", **extra}
    par.write_text("".join(f"{k}: {v}\n" for k, v in keys.items()))
    return par


def subset() -> Path:
    """convertf's PACKEDANCESTRYMAP subset of AADR v66.p1 to POPS (made once)."""
    out = CACHE / "bench16"
    if not Path(f"{out}.geno").exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        pops = CACHE / "bench16.pops"
        pops.write_text("".join(f"{p}\n" for p in POPS))
        par = _convertf_par(AADR, out, "PACKEDANCESTRYMAP", poplistname=str(pops))
        # convertf writes a profile (gmon.out) where it runs
        subprocess.run([CONVERTF, "-p", str(par)], cwd=CACHE, check=True, stdout=subprocess.DEVNULL)
    return out


@dataclass
class Case:
    name: str
    analysis: str
    reference: str
    popmix: Callable[[Path, Path, int], list[str]]  # (input, outdir, threads) -> argv
    ref: Callable[[Path, Path], list[str]]  # (input, outdir) -> argv
    agree: Callable[[Path, Path], dict]  # (popmix outdir, reference outdir) -> agreement
    prepare: Callable[[], Path] = subset


# --- convertf ---------------------------------------------------------------

def _eigenstrat_rows(prefix: Path) -> dict[str, tuple[str, str, str]]:
    snps = [line.split() for line in Path(f"{prefix}.snp").read_text().splitlines()]
    with open(f"{prefix}.geno") as g:
        return {s[0]: (s[4], s[5], line.rstrip("\n")) for s, line in zip(snps, g)}


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


# --- f2 extraction ----------------------------------------------------------

_R_F2 = """
suppressMessages(library(admixtools))
a <- commandArgs(TRUE)
pops <- strsplit(a[2], ",")[[1]]
f2 <- f2_from_geno(a[1], pops = pops, verbose = FALSE)
writeBin(as.vector(f2), file.path(a[3], "f2.bin"))
jsonlite::write_json(list(dim = dim(f2), pops = dimnames(f2)[[1]]), file.path(a[3], "f2.json"),
                     auto_unbox = TRUE)
"""


def _r(script: str, outdir: Path, *args: str) -> list[str]:
    path = outdir / "run.R"
    path.write_text(script)
    return [RSCRIPT, str(path), *args, str(outdir)]


def _agree_f2(ours: Path, theirs: Path) -> dict:
    import popmix as pm

    meta = json.loads((theirs / "f2.json").read_text())
    n, _, nb = meta["dim"]
    ref = np.fromfile(theirs / "f2.bin", dtype="<f8").reshape((nb, n, n)).transpose(1, 2, 0)
    f2 = pm.read_f2(str(ours / "f2"))
    worst_abs = worst_rel = 0.0
    for i, a in enumerate(meta["pops"]):
        for j, b in enumerate(meta["pops"]):
            if i < j:
                got, want = np.asarray(f2.f2_per_block(a, b)), ref[i, j]
                d = np.abs(got - want)
                worst_abs = max(worst_abs, float(np.nanmax(d)))
                worst_rel = max(worst_rel, float(np.nanmax(d / np.maximum(np.abs(want), 1e-300))))
    ok = list(f2.block_lengths) and len(f2.block_lengths) == nb
    return {"blocks": nb, "max_abs": worst_abs, "max_rel": worst_rel,
            "pass": bool(ok) and (worst_rel <= 1e-10 or worst_abs <= 1e-15)}


F2_CASE = Case(
    name="extract_f2", analysis=f"f2 blocks, {len(POPS)} populations", reference="admixtools2 2.0.10 f2_from_geno",
    popmix=lambda i, o, t: ["popmix", "--threads", str(t), "extract-f2", "-i", str(i), "-o", str(o / "f2"),
                            "--populations", ",".join(POPS)],
    ref=lambda i, o: _r(_R_F2, o, str(i), ",".join(POPS)),
    agree=_agree_f2,
)


# --- qpAdm from genotypes ---------------------------------------------------

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


def _popmix_qpadm(i: Path, o: Path, t: int) -> list[str]:
    q = QPADM
    return ["sh", "-c", " ".join(["popmix", "--threads", str(t), "qpadm", "-i", str(i), "--target", q["target"],
                                  "--left", ",".join(q["left"]), "--right", ",".join(q["right"]),
                                  "--format", "json", ">", str(o / "qpadm.json")])]


QPADM_CASE = Case(
    name="qpadm", analysis="qpAdm from genotypes (2 sources, 7 right)", reference="admixtools2 2.0.10 qpadm",
    popmix=_popmix_qpadm,
    ref=lambda i, o: _r(_R_QPADM, o, str(i), QPADM["target"], ",".join(QPADM["left"]), ",".join(QPADM["right"])),
    agree=_agree_qpadm,
)

# --- PCA ------------------------------------------------------------------------

SMARTPCA = os.environ.get("SMARTPCA_HIPREC",
                          str(HOME / ".local/share/popmix-references/eigensoft-7.2.1-hiprec/smartpca"))
PLINK2 = os.environ.get("POPMIX_PLINK2", "plink2")
PRESENT_DAY = ["Mbuti", "Yoruba", "French", "Han", "Papuan", "Karitiana", "Basque", "Sardinian"]


def present_day() -> Path:
    """The subset with every sample outside PRESENT_DAY labelled Ignore, so
    neither tool uses or projects it."""
    src = subset()
    out = CACHE / "bench16_present"
    if not Path(f"{out}.ind").exists():
        for ext in (".geno", ".snp"):
            Path(f"{out}{ext}").symlink_to(f"{src}{ext}")
        keep = set(PRESENT_DAY)
        lines = [line.split() for line in Path(f"{src}.ind").read_text().splitlines()]
        Path(f"{out}.ind").write_text("".join(f"{i} {s} {p if p in keep else 'Ignore'}\n" for i, s, p in lines))
    return out


def _smartpca(i: Path, o: Path) -> list[str]:
    (o / "pops.txt").write_text("".join(f"{p}\n" for p in PRESENT_DAY))
    par = o / "smartpca.par"
    par.write_text("".join(f"{k}: {v}\n" for k, v in {
        "genotypename": f"{i}.geno", "snpname": f"{i}.snp", "indivname": f"{i}.ind",
        "evecoutname": o / "out.evec", "evaloutname": o / "out.eval", "poplistname": o / "pops.txt",
        "numoutevec": 10, "numoutlieriter": 0}.items()))
    return [SMARTPCA, "-p", str(par)]


def _agree_pca(ours: Path, theirs: Path) -> dict:
    got = np.array([float(x) for x in (ours / "out.eval").read_text().split()][:10])
    want = np.array([float(x) for x in (theirs / "out.eval").read_text().split()][:10])
    rel = float(np.max(np.abs(got - want) / np.abs(want)))
    return {"eigenvalues": len(want), "max_rel": rel, "pass": rel <= 1e-9}


PCA_CASE = Case(
    name="pca", analysis=f"PCA, {len(PRESENT_DAY)} present-day populations, 10 PCs",
    reference="EIGENSOFT 7.2.1 smartpca",
    popmix=lambda i, o, t: ["popmix", "--threads", str(t), "pca", "-i", str(i), "--pca-pops", ",".join(PRESENT_DAY),
                            "--components", "10", "--outlier-iterations", "0", "-o", str(o / "out")],
    ref=_smartpca, agree=_agree_pca, prepare=present_day,
)


# --- KING-robust kinship -----------------------------------------------------

def present_day_bed() -> Path:
    """convertf's PLINK .bed of the PRESENT_DAY samples."""
    src = subset()
    out = CACHE / "bench16_present_bed"
    if not Path(f"{out}.bed").exists():
        pops = CACHE / "present.pops"
        pops.write_text("".join(f"{p}\n" for p in PRESENT_DAY))
        par = CACHE / "present_bed.par"
        par.write_text(f"genotypename: {src}.geno\nsnpname: {src}.snp\nindivname: {src}.ind\n"
                       f"outputformat: PACKEDPED\ngenotypeoutname: {out}.bed\nsnpoutname: {out}.bim\n"
                       f"indivoutname: {out}.fam\npoplistname: {pops}\n")
        subprocess.run([CONVERTF, "-p", str(par)], cwd=CACHE, check=True, stdout=subprocess.DEVNULL)
    return out


def _agree_king(ours: Path, theirs: Path) -> dict:
    ids = [line.split("\t")[-1] for line in (theirs / "king.king.id").read_text().splitlines()[1:]]
    tri = np.fromfile(theirs / "king.king.bin", dtype=np.float64)
    want, k = {}, 0
    for a in range(1, len(ids)):
        for b in range(a):
            want[frozenset((ids[a], ids[b]))] = tri[k]
            k += 1
    rows = [line.split("\t") for line in (ours / "kin.tsv").read_text().splitlines()[1:]]
    def diff(got: float, ref: float) -> float:
        # a sample with no heterozygous call gives -inf in both
        return 0.0 if got == ref else abs(got - ref)

    worst = max(diff(float(r[-1]), want[frozenset((r[0], r[1]))]) for r in rows)
    return {"pairs": len(rows), "max_abs": worst, "pass": len(rows) == len(want) and worst <= 1e-10}


KING_CASE = Case(
    name="kinship", analysis=f"KING-robust kinship, every pair of {len(PRESENT_DAY)} populations' samples",
    reference="PLINK 2 --make-king-table (16 threads)",
    popmix=lambda i, o, t: ["popmix", "--threads", str(t), "kinship", "-i", str(i), "--method", "king",
                            "-o", str(o / "kin.tsv")],
    ref=lambda i, o: [PLINK2, "--bfile", str(i), "--make-king-table", "--make-king", "triangle", "bin",
                      "--threads", "16", "--out", str(o / "king")],
    agree=_agree_king, prepare=present_day_bed,
)

# --- qpGraph -----------------------------------------------------------------

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
r <- qpgraph(f2, edges, numstart = 10, seed = 1, verbose = FALSE)
jsonlite::write_json(list(score = r$score), file.path(a[4], "qpgraph.json"), digits = NA, auto_unbox = TRUE)
"""


def _graph_file(o: Path) -> Path:
    path = o / "graph.txt"
    path.write_text("".join(f"{a}\t{b}\n" for a, b in GRAPH))
    return path


def _agree_qpgraph(ours: Path, theirs: Path) -> dict:
    got = json.loads((ours / "qpgraph.json").read_text())["score"]
    want = json.loads((theirs / "qpgraph.json").read_text())["score"]
    rel = abs(got - want) / abs(want)
    # an optimized score: R's L-BFGS-B stops at factr * epsilon
    return {"score": want, "score_rel": rel, "pass": rel <= 1e-6}


QPGRAPH_CASE = Case(
    name="qpgraph", analysis="qpGraph, 6 populations, 1 admixture, 10 starts",
    reference="admixtools2 2.0.10 qpgraph",
    popmix=lambda i, o, t: ["sh", "-c", f"popmix --threads {t} qpgraph -i {i} --graph {_graph_file(o)} "
                                        f"--numstart 10 --seed 1 --format json > {o / 'qpgraph.json'}"],
    ref=lambda i, o: _r(_R_QPGRAPH, o, str(i), str(_graph_file(o)), ",".join(GRAPH_POPS)),
    agree=_agree_qpgraph,
)


# --- qpDstat from genotypes --------------------------------------------------

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

CASES = {c.name: c for c in (CONVERTF_CASE, F2_CASE, QPADM_CASE, PCA_CASE, KING_CASE, QPGRAPH_CASE,
                             QPDSTAT_CASE)}
