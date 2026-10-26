"""Benchmark cases: popmix and a reference tool on identical inputs.

Every case runs the popmix command line (one process, as a user would) and
the reference tool's own program, each writing its result to its own
directory; `agree` then compares the two results at the CLAUDE.md
tolerances. Inputs are made once, by the reference tools, in the cache.
"""

from __future__ import annotations

import json
import os
import re
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


def _f2_mismatches(got: np.ndarray, want: np.ndarray) -> int:
    """Entries off by more than CLAUDE.md's f-statistic tolerance, atol 1e-15 + rtol 1e-10.

    NaN must meet NaN.
    """
    nan = np.isnan(got) | np.isnan(want)
    off = np.abs(got - want) > 1e-15 + 1e-10 * np.abs(want)
    return int(np.sum(np.where(nan, np.isnan(got) != np.isnan(want), off)))


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
r <- qpgraph(f2, edges, numstart = 10, seed = as.integer(a[4]), verbose = FALSE)
jsonlite::write_json(list(score = r$score), file.path(a[5], "qpgraph.json"), digits = NA, auto_unbox = TRUE)
"""


def _graph_file(o: Path) -> Path:
    path = o / "graph.txt"
    path.write_text("".join(f"{a}\t{b}\n" for a, b in GRAPH))
    return path


def _score(name: str) -> Callable[[Path], float]:
    return lambda d: json.loads((d / name).read_text())["score"]


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

# --- ADMIXTURE ------------------------------------------------------------------

ADMIXTURE = os.environ.get("POPMIX_ADMIXTURE", "/usr/local/bin/admixture")
ADMIXTURE_POPS = ["Yoruba", "French", "Han"]
ADMIXTURE_K = 3
# as tests/golden/admixture: tight convergence, so each tool stops close to
# the maximum it climbs to
ADMIXTURE_TOL = 1e-10
SPREAD_SEEDS = (2, 3)


def admixture_bed() -> Path:
    """PLINK 2's LD-pruned autosomal .bed of the ADMIXTURE_POPS samples
    (SNPs called in at least 90% of them, MAF >= 1%; r2 < 0.1 in 200-SNP
    windows), from the present-day .bed."""
    src = present_day_bed()
    out = CACHE / "bench16_admixture"
    if not Path(f"{out}.bed").exists():
        pop = {line.split()[0]: line.split()[2] for line in Path(f"{subset()}.ind").read_text().splitlines()}
        keep = CACHE / "admixture.keep"
        keep.write_text("".join(f"{f} {i}\n" for f, i, *_ in (line.split() for line in
                                                                 Path(f"{src}.fam").read_text().splitlines())
                                if pop[i] in ADMIXTURE_POPS))
        common = [PLINK2, "--bfile", str(src), "--keep", str(keep), "--threads", "4"]
        subprocess.run([*common, "--autosome", "--geno", "0.1", "--maf", "0.01", "--indep-pairwise", "200", "50",
                        "0.1", "--out", str(CACHE / "admixture_prune")], check=True, stdout=subprocess.DEVNULL)
        subprocess.run([*common, "--extract", str(CACHE / "admixture_prune.prune.in"), "--make-bed", "--out",
                        str(out)], check=True, stdout=subprocess.DEVNULL)
    return out


def _popmix_admixture(i: Path, o: Path, t: int, seed: int = 1) -> list[str]:
    return ["popmix", "--threads", str(t), "admixture", "-i", str(i), "-k", str(ADMIXTURE_K), "--seed", str(seed),
            "--tolerance", str(ADMIXTURE_TOL), "-o", str(o / "out")]


def _admixture(i: Path, o: Path, seed: int = 1) -> list[str]:
    # writes <input name>.K.Q and .P where it runs
    return [ADMIXTURE, f"--seed={seed}", f"-C={ADMIXTURE_TOL}", f"{i}.bed", str(ADMIXTURE_K)]


def _loglik(log: Path) -> float:
    """The last log-likelihood either tool printed (both print 'Loglikelihood: x')."""
    return float(re.findall(r"^Loglikelihood: (\S+)", log.read_text(), re.M)[-1])


def _seed_runs(theirs: Path, tool: str, cmd: Callable[[Path, int], list[str]],
               value: Callable[[Path], float]) -> dict[int, float]:
    """The result of each SPREAD_SEEDS run of one tool (run once per
    benchmark, kept with the reference output): the tool's run-to-run spread."""
    cached = theirs / f"seeds-{tool}.json"
    if not cached.exists():
        out = {}
        for seed in SPREAD_SEEDS:
            d = theirs / f"seed{seed}-{tool}"
            d.mkdir(exist_ok=True)
            run = subprocess.run(cmd(d, seed), cwd=d, check=True, capture_output=True, text=True)
            (d / "run.log").write_text(run.stdout)
            out[seed] = value(d)
        cached.write_text(json.dumps(out))
    return {int(k): v for k, v in json.loads(cached.read_text()).items()}


def _seed_logliks(theirs: Path, tool: str, cmd: Callable[[Path, int], list[str]]) -> dict[int, float]:
    return _seed_runs(theirs, tool, cmd, lambda d: _loglik(d / "run.log"))


def _agree_admixture(ours: Path, theirs: Path) -> dict:
    bed = admixture_bed()
    got, want = _loglik(ours / "run.log"), _loglik(theirs / "run.log")
    mine = [got, *_seed_logliks(theirs, "popmix", lambda d, s: _popmix_admixture(bed, d, 4, s)).values()]
    ref = [want, *_seed_logliks(theirs, "admixture", lambda d, s: _admixture(bed, d, s)).values()]
    # as the golden test: popmix's maximum is at least ADMIXTURE's
    return {"loglik": got, "reference_loglik": want, "loglik_diff": got - want,
            "seeds": 1 + len(SPREAD_SEEDS), "popmix_spread": max(mine) - min(mine),
            "reference_spread": max(ref) - min(ref), "pass": got >= want - 1e-6}


ADMIXTURE_CASE = Case(
    name="admixture", analysis=f"ADMIXTURE K={ADMIXTURE_K}, {len(ADMIXTURE_POPS)} present-day populations, "
                               "LD-pruned, seed 1",
    reference="ADMIXTURE 1.3.0 (1 thread)",
    popmix=_popmix_admixture, ref=_admixture, agree=_agree_admixture, prepare=admixture_bed,
)


# --- ROLLOFF -------------------------------------------------------------------

CLASSIC = Path(os.environ.get("POPMIX_ADMIXTOOLS7", HOME / ".local/share/popmix-references/admixtools-7.0.2"))
# two chromosomes (the jackknife needs two), to 5 cM: rolloff takes about 7
# minutes for chromosome 22 alone, and over half an hour for 20-22 to 10 cM
ROLLOFF = {"admixed": "French", "refs": ["Sardinian", "Russia_Samara_EBA_Yamnaya"], "binsize": 0.001,
           "maxdis": 0.05, "chromosomes": ["21", "22"]}


def _rolloff_badsnps_path() -> Path:
    return CACHE / f"rolloff_chr{'_'.join(ROLLOFF['chromosomes'])}.badsnps"


def rolloff_badsnps() -> Path:
    """The subset's SNPs off ROLLOFF's chromosomes (rolloff's badsnpname, as
    generate_rolloff.py). Returns the subset."""
    src = subset()
    bad = _rolloff_badsnps_path()
    if not bad.exists():
        keep = set(ROLLOFF["chromosomes"])
        bad.write_text("".join(f"{s[0]}\n" for s in (line.split() for line in
                                                      Path(f"{src}.snp").read_text().splitlines())
                               if s[1] not in keep))
    return src


def _rolloff(i: Path, o: Path) -> list[str]:
    # tests/golden/rolloff/generate_rolloff.py's parameters
    (o / "refs").write_text("".join(f"{p}\n" for p in ROLLOFF["refs"]))
    par = o / "rolloff.par"
    par.write_text("".join(f"{k}: {v}\n" for k, v in {
        "genotypename": f"{i}.geno", "snpname": f"{i}.snp", "indivname": f"{i}.ind", "poplistname": o / "refs",
        "admixpop": ROLLOFF["admixed"], "output": o / "rolloff.out", "binsize": ROLLOFF["binsize"],
        "maxdis": ROLLOFF["maxdis"], "runmode": 1, "chithresh": 6.0, "zdipcorrmode": "YES", "checkmap": "YES",
        "seed": 1, "numchrom": 22, "jackknife": "YES", "badsnpname": _rolloff_badsnps_path()}.items()))
    return [str(CLASSIC / "rolloff-hiprec"), "-p", str(par)]


def _agree_rolloff(ours: Path, theirs: Path) -> dict:
    want = np.array([[float(x) for x in line.split()] for line in (theirs / "rolloff.out").read_text().splitlines()
                     if line.strip() and not line.lstrip().startswith("#")])
    import pandas as pd

    got = pd.read_csv(ours / "out.bins.tsv", sep="\t")
    # tests/golden/test_rolloff.py: curves rtol 1e-9, pair counts exact
    cols = {"zscore": 1, "regression": 2, "corr": 3}
    rel = {k: float(np.max(np.abs(got[k].to_numpy() - want[:, j]) / np.maximum(np.abs(want[:, j]), 1e-12)))
           for k, j in cols.items()} if len(got) == len(want) else {k: float("inf") for k in cols}
    counts = (len(got) == len(want) and bool((got["count"].to_numpy() == want[:, 4]).all())
              and np.allclose(got["d_cm"].to_numpy(), want[:, 0], rtol=1e-12, atol=0))
    return {"bins": len(want), **{f"{k}_rel": v for k, v in rel.items()}, "bins_equal": counts,
            "pass": counts and max(rel.values()) <= 1e-9}


ROLLOFF_CASE = Case(
    name="rolloff", analysis=f"ROLLOFF, {ROLLOFF['admixed']} = {' + '.join(ROLLOFF['refs'])}, "
                             f"to {ROLLOFF['maxdis'] * 100:g} cM, chromosomes {ROLLOFF['chromosomes'][0]}-"
                             f"{ROLLOFF['chromosomes'][-1]}, jackknife over them",
    reference="ADMIXTOOLS 7.0.2 rolloff (full-precision build)",
    popmix=lambda i, o, t: ["popmix", "--threads", str(t), "rolloff", "-i", str(i), "--admixed", ROLLOFF["admixed"],
                            "--refs", ",".join(ROLLOFF["refs"]), "--binsize", str(ROLLOFF["binsize"]),
                            "--maxdis", str(ROLLOFF["maxdis"]), "--chromosomes", ",".join(ROLLOFF["chromosomes"]),
                            "-o", str(o / "out")],
    ref=_rolloff, agree=_agree_rolloff, prepare=rolloff_badsnps,
)


# --- DATES ---------------------------------------------------------------------

DATES = Path(os.environ.get("POPMIX_DATES", HOME / "DATES"))
# Russia_Samara_EBA_Yamnaya mixes pseudo-haploid samples and a diploid genome
# (I6729.DG): --compat dates-v753 counts every sample diploid, as DATES does
# (DEVIATIONS.md, DATES-0004)
DATES_POPS = {"test": "French", "source1": "Sardinian", "source2": "Russia_Samara_EBA_Yamnaya"}
DATES_BINSIZE, DATES_MAXDIS, DATES_QBIN, DATES_LOVAL = 0.001, 0.5, 10, 0.45
# tests/golden/test_dates_v753.py: printed precision, plus dates_expfit's
# optimizer tolerance (gslprecision 1e-6 in r) for fitted generations
COV_TOL, CORR_TOL, GEN_TOL = 5e-10 + 1e-15, 5e-7 + 1e-12, 5e-4 + 1e-9
FIT_TOL = GEN_TOL + 1e-6 / DATES_BINSIZE


REPO = Path(__file__).resolve().parents[1]
DATES_BUILD = CACHE / "dates-v753-patched"


def dates_tools() -> Path:
    """DATES v753 built with tests/golden/dates/dates_numbins_fix.patch, as the
    goldens use, and its helper scripts on one PATH (DATES_BUILD/tools).

    The stock bin/dates reads an uninitialized buffer size and aborts at
    random ("Unable to allocate -1817970720 unit(s)"), in mesh mode too, and
    its dowtjack needs a libgsl this system no longer has. The build is from
    the DATES checkout's committed source, in a clean environment (conda's
    compiler flags break the link). Returns the subset."""
    tools = DATES_BUILD / "tools"
    if not (tools / "dowtjack").exists():
        src = DATES_BUILD / "src"
        DATES_BUILD.mkdir(parents=True, exist_ok=True)
        archive = subprocess.run(["git", "-C", str(DATES), "archive", "HEAD", "src"], check=True,
                                 capture_output=True).stdout
        subprocess.run(["tar", "-x", "-C", str(DATES_BUILD)], input=archive, check=True)
        # the archive carries prebuilt objects and libraries: build from source
        subprocess.run(["sh", "-c", "rm -f *.o nicksrc/*.o nicksrc/libnick.a nicklib dates dowtjack"], cwd=src,
                       check=True)
        subprocess.run(["patch", "-p1", "-i", str(REPO / "tests/golden/dates/dates_numbins_fix.patch")], cwd=src,
                       check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["make", "dates", "dowtjack"], cwd=src, check=True, capture_output=True,
                       env={"PATH": "/usr/local/bin:/usr/bin:/bin", "HOME": str(HOME)})
        tools.mkdir(exist_ok=True)
        for prog in [*(DATES / "bin").iterdir(), src / "dates", src / "dowtjack"]:
            link = tools / prog.name
            if link.is_symlink():
                link.unlink()
            link.symlink_to(prog)
    return subset()


def _dates(i: Path, o: Path) -> list[str]:
    p = DATES_POPS
    (o / "admixlist").write_text(f"{p['source1']} {p['source2']} {p['test']} {p['test']}\n")
    par = o / "dates.par"
    # tests/golden/dates/generate_dates_reference.sh's parameters
    par.write_text("".join(f"{k}: {v}\n" for k, v in {
        "genotypename": f"{i}.geno", "snpname": f"{i}.snp", "indivname": f"{i}.ind", "admixlist": "admixlist",
        "binsize": DATES_BINSIZE, "maxdis": DATES_MAXDIS, "seed": 77, "runmode": 1, "chithresh": 0.0,
        "mincount": 4, "qbin": DATES_QBIN, "runfit": "YES", "afffit": "YES", "lovalfit": DATES_LOVAL}.items()))
    tools = DATES_BUILD / "tools"
    # dates runs its helpers (dates_jackknife, dowtjack, ...) by name
    return ["env", f"PATH={tools}:{os.environ['PATH']}", str(tools / "dates"), "-p", par.name]


def _agree_dates(ours: Path, theirs: Path) -> dict:
    test = DATES_POPS["test"]
    res = theirs / test
    got = json.loads((ours / "dates.json").read_text())
    rows = [line.split() for line in (res / f"{test}.out").read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    # the JSON holds the first 20 bins
    ld = [(abs(b["n_pairs"] - int(r[4])), abs(b["weighted_ld"] - float(r[1])), abs(b["regression_coef"] - float(r[2])),
           abs(b["correlation"] - float(r[3]))) for b, r in zip(got["ld_bins"], rows)]  # noqa: B905
    full = float(re.search(r"mean \(generations\):\s+(\S+)", (res / f"expfit_{test}.flog").read_text()).group(1))
    blocks = sum(float(line.split()[1]) > 0 for line in (res / f"{test}.jin").read_text().splitlines())
    jk_est, jk_se = (float(x) for x in (res / f"{test}.jout").read_text().split())
    jk_tol = (blocks + 2) * FIT_TOL + GEN_TOL
    out = {"bins_compared": len(ld), "pairs_equal": all(d[0] == 0 for d in ld),
           "max_cov_abs": max(d[1] for d in ld), "max_corr_abs": max(max(d[2], d[3]) for d in ld),
           "fit_abs": abs(got["fit"]["decay_rate"] - full), "generations_abs": abs(got["generations"] - jk_est),
           "se_abs": abs(got["generations_se"] - jk_se)}
    out["pass"] = (out["pairs_equal"] and out["max_cov_abs"] <= COV_TOL and out["max_corr_abs"] <= CORR_TOL
                   and out["fit_abs"] <= FIT_TOL and out["generations_abs"] <= jk_tol and out["se_abs"] <= jk_tol)
    return out


def _popmix_dates(i: Path, o: Path, t: int) -> list[str]:
    p = DATES_POPS
    return ["sh", "-c", f"popmix --threads {t} dates -i {i} --test-pop {p['test']} --source1 {p['source1']} "
                        f"--source2 {p['source2']} --bin-width {DATES_BINSIZE * 100:g} "
                        f"--max-distance {DATES_MAXDIS * 100:g} --min-distance {DATES_LOVAL} --qbin {DATES_QBIN} "
                        f"--compat dates-v753 --format json > {o / 'dates.json'}"]


DATES_CASE = Case(
    name="dates", analysis=f"DATES, {DATES_POPS['test']} = {DATES_POPS['source1']} + {DATES_POPS['source2']}, "
                           f"FFT mesh (qbin {DATES_QBIN}), to {DATES_MAXDIS * 100:g} cM, jackknife",
    reference="DATES v753 (numbins fix, as the goldens)",
    popmix=_popmix_dates, ref=_dates, agree=_agree_dates, prepare=dates_tools,
)


# --- READ kinship ----------------------------------------------------------------

READV2 = Path(os.environ.get("POPMIX_READV2", HOME / ".local/share/popmix-references/readv2"))
READ_POP = "Spain_EN"


def read_bed() -> Path:
    """convertf's PLINK .bed of the READ_POP samples."""
    src = subset()
    out = CACHE / "bench16_read_bed"
    if not Path(f"{out}.bed").exists():
        pops = CACHE / "read.pops"
        pops.write_text(f"{READ_POP}\n")
        par = CACHE / "read_bed.par"
        par.write_text(f"genotypename: {src}.geno\nsnpname: {src}.snp\nindivname: {src}.ind\n"
                       f"outputformat: PACKEDPED\ngenotypeoutname: {out}.bed\nsnpoutname: {out}.bim\n"
                       f"indivoutname: {out}.fam\npoplistname: {pops}\n")
        subprocess.run([CONVERTF, "-p", str(par)], cwd=CACHE, check=True, stdout=subprocess.DEVNULL)
    return out


def _readv2(i: Path, o: Path) -> list[str]:
    # READv2 writes a filtered copy of its input next to it: give it links
    for ext in (".bed", ".bim", ".fam"):
        link = o / f"data{ext}"
        if not link.exists():
            link.symlink_to(f"{i}{ext}")
    return [str(READV2 / "venv/bin/python"), str(READV2 / "READ2.py"), "-i", "data"]


# tests/golden/test_read_readv2.py: popmix's column, READv2's
READ_COLUMNS = {"normalized_p0": "P0_mean", "p0": "Nonnormalized_P0", "p0_se": "Nonnormalized_P0_serr",
                "first_degree_windows_outside": "Perc_Win_1stdeg_P0", "effective_snps": "NSNPsXNorm",
                "kinship": "KinshipCoefficient", "z_upper": "Zup", "z_lower": "Zdown"}


def _agree_read(ours: Path, theirs: Path) -> dict:
    import pandas as pd

    want = pd.read_csv(theirs / "Read_Results.tsv", sep="\t")
    got = pd.read_csv(ours / "kin.tsv", sep="\t", keep_default_na=False, na_values=["", "nan", "NaN"])
    got.index = got["id1"] + "," + got["id2"]
    if set(got.index) != set(want["PairIndividuals"]):
        return {"pairs": len(want), "pass": False}
    got = got.loc[want["PairIndividuals"]]
    # relationship and first-degree type exact (as in the golden test)
    classes = (got["relationship"].tolist() == want["Rel"].tolist()
               and got["first_degree_type"].fillna("").tolist() == want["1st_Type"].fillna("").tolist()
               and got["overlap_snps"].tolist() == want["OverlapNSNPs"].tolist())
    worst = 0.0
    for ours_col, theirs_col in READ_COLUMNS.items():
        a, b = got[ours_col].to_numpy(dtype=float), want[theirs_col].to_numpy(dtype=float)
        both = np.isnan(a) & np.isnan(b)
        worst = max(worst, float(np.max(np.where(both, 0.0, np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))))
    return {"pairs": len(want), "classes_equal": classes, "max_rel": worst, "pass": classes and worst <= 1e-10}


READ_CASE = Case(
    name="read_kinship", analysis=f"READ kinship, every pair of {READ_POP} samples", reference="READv2 (73739be)",
    popmix=lambda i, o, t: ["popmix", "--threads", str(t), "kinship", "-i", str(i), "--method", "read",
                            "--compat", "readv2-2.01", "-o", str(o / "kin.tsv")],
    ref=_readv2, agree=_agree_read, prepare=read_bed,
)


# --- classic ADMIXTOOLS 7.0.2 -------------------------------------------------
# Classic ADMIXTOOLS prints a few decimals, so agreement is to half a unit in
# the last printed digit. With inbreed: NO (its default) it treats every
# sample as diploid, so popmix runs with --pseudohaploid diploid.

CLASSIC_BIN = CLASSIC / "bin"
PRINTED6, PRINTED3 = 5e-7 + 1e-12, 5e-4 + 1e-9


def _classic_par(i: Path, o: Path, name: str, **keys) -> Path:
    par = o / f"{name}.par"
    par.write_text("".join(f"{k}: {v}\n" for k, v in {"genotypename": f"{i}.geno", "snpname": f"{i}.snp",
                                                       "indivname": f"{i}.ind", **keys}.items()))
    return par


def _classic(prog: str, par: Path, *args: str) -> list[str]:
    # the classic programs print their results to stdout
    return ["sh", "-c", f"{CLASSIC_BIN / prog} -p {par} {' '.join(args)} > {par.with_suffix('.out')}"]


def _results(out: Path) -> list[list[str]]:
    return [line.split()[1:] for line in out.read_text().splitlines() if line.split()[:1] == ["result:"]]


def _quads() -> list[tuple[str, ...]]:
    return [(a, b, c, d) for a in DSTAT["pop1"] for b in DSTAT["pop2"] for c in DSTAT["pop3"] for d in DSTAT["pop4"]]


def _classic_qpdstat(i: Path, o: Path) -> list[str]:
    (o / "quads").write_text("".join(" ".join(q) + "\n" for q in _quads()))
    return _classic("qpDstat", _classic_par(i, o, "qpdstat", popfilename=o / "quads", f4mode="YES"))


def _agree_classic_qpdstat(ours: Path, theirs: Path) -> dict:
    got = {tuple(r[k] for k in ("pop1", "pop2", "pop3", "pop4")): r
           for r in json.loads((ours / "qpdstat.json").read_text())}
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
                                       for c in QP3POP["pop3"]))
    return _classic("qp3Pop", _classic_par(i, o, "qp3pop", popfilename=o / "triples", outgroupmode="YES"))


def _agree_classic_qp3pop(ours: Path, theirs: Path) -> dict:
    got = {(r["pop1"], r["pop2"], r["pop3"]): r for r in json.loads((ours / "qp3pop.json").read_text())}
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
    (o / "left").write_text("".join(f"{p}\n" for p in [q["target"], *q["left"]]))
    (o / "right").write_text("".join(f"{p}\n" for p in q["right"]))
    return _classic("qpAdm", _classic_par(i, o, "qpadm", popleft=o / "left", popright=o / "right"))


def _agree_classic_qpadm(ours: Path, theirs: Path) -> dict:
    got = json.loads((ours / "qpadm.json").read_text())
    text = (theirs / "qpadm.out").read_text()
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
    path.write_text("\n".join(lines) + "\n")
    return path


def _classic_qpgraph(i: Path, o: Path) -> list[str]:
    # tests/golden/admixtools_classic (removed in 9ed980b) ran qpGraph so
    par = _classic_par(i, o, "qpgraph", outpop="NULL", useallsnps="NO", blgsize=0.05, forcezmode="YES",
                       lsqmode="NO", diag=".0001", bigiter=6, hires="YES", lambdascale=1)
    return _classic("qpGraph", par, "-g", str(_classic_graph(o)), "-o", str(o / "qpgraph.ggg"))


def _agree_classic_qpgraph(ours: Path, theirs: Path) -> dict:
    got = json.loads((ours / "qpgraph.json").read_text())["score"]
    want = float(re.search(r"final score:\s+(\S+)", (theirs / "qpgraph.out").read_text()).group(1))
    # printed to 3 decimals; both optimize, R's L-BFGS-B to about 1e-8 of the score
    return {"score": want, "score_abs": abs(got - want), "pass": abs(got - want) <= PRINTED3 + 1e-6 * abs(want)}


CLASSIC_QPGRAPH_CASE = Case(
    name="classic_qpgraph", analysis="qpGraph, 6 populations, 1 admixture", reference="ADMIXTOOLS 7.0.2 qpGraph",
    popmix=lambda i, o, t: ["sh", "-c", f"popmix --threads {t} qpgraph -i {i} --graph {_graph_file(o)} "
                                        f"--pseudohaploid diploid --numstart 10 --seed 1 --format json "
                                        f"> {o / 'qpgraph.json'}"],
    ref=_classic_qpgraph, agree=_agree_classic_qpgraph,
)


# --- admixtools2: qpAdm rotation, find_graphs, resampling ----------------------
# popmix has no command for rotation or resampling: its side is one Python
# process with the popmix API, as the CLI would be.

def _py(script: str, outdir: Path, threads: int, *args: str) -> list[str]:
    path = outdir / "run.py"
    path.write_text(script)
    return ["python", str(path), str(threads), *args, str(outdir)]


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
    got = json.loads((ours / "find_graphs.json").read_text())
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

CASES = {c.name: c for c in (CONVERTF_CASE, F2_CASE, QPADM_CASE, PCA_CASE, KING_CASE, QPGRAPH_CASE,
                             QPDSTAT_CASE, ADMIXTURE_CASE, ROLLOFF_CASE, DATES_CASE, READ_CASE,
                             CLASSIC_QPDSTAT_CASE, CLASSIC_QP3POP_CASE, CLASSIC_QPADM_CASE, CLASSIC_QPGRAPH_CASE,
                             ROTATE_CASE, FIND_GRAPHS_CASE, RESAMPLE_CASE)}
