"""Benchmark cases: smartpca, KING, ADMIXTURE and READ."""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path
from typing import Callable

import numpy as np

from ._common import (
    CACHE,
    CONVERTF,
    HOME,
    PLINK2,
    PRESENT_DAY,
    SMARTPCA,
    Case,
    present_day,
    present_day_bed,
    subset,
)


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
