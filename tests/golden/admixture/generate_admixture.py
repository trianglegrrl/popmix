#!/usr/bin/env python3
"""Generate the ADMIXTURE 1.3.0 reference for popmix's admixture.

Usage: python tests/golden/admixture/generate_admixture.py [config...]

Each config's samples are written as PLINK files by EIGENSOFT convertf
(restricted to some chromosomes by PLINK 2 where the config says, dropping
SNPs with no genotype at all, which ADMIXTURE refuses), and
ADMIXTURE 1.3.0 is run with --seed=1 and -C=1e-10 (tight convergence, so
its estimates are close to the maximum it converges to) in one of three
modes:

  - unsupervised, with 5-fold cross-validation (--cv=5)
  - projection (-P): Q for new samples with P fixed to that of an
    unsupervised run on other samples (the run's .P as <prefix>.K.P.in)
  - supervised (--supervised, a .pop file labeling the reference samples)

and tests/golden/admixture/v1.3.0/<config>/admixture.json records the sample
IDs (.fam order), K, the final log-likelihood and the iterations it reports,
the CV error, Q and P as written (6 decimals), and for projection the P it
was given.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AADR = "/references/aadr/v62.0_1240k_public"
SMALL = str(REPO / "tests/golden/dates/family_chr20_22")
ADMIXTURE = os.environ.get("ADMIXTURE", "admixture")
PLINK2 = os.environ.get("PLINK2", "plink2")
CONVERTF = os.environ.get("CONVERTF", "/usr/bin/convertf")
FLAGS = ["--seed=1", "-C=1e-10"]

KGP = ["CEU.DG", "YRI.DG", "CHB.DG"]
CONFIGS = {
    "family_k2": {"prefix": SMALL, "pops": ["French", "Yoruba", "Simulation", "Simulation.1"], "k": 2,
                      "mode": "unsupervised"},
    "family_k3": {"prefix": SMALL, "pops": ["French", "Yoruba", "Simulation", "Simulation.1"], "k": 3,
                      "mode": "unsupervised"},
    "family_projection": {"prefix": SMALL, "pops": ["French", "Yoruba"], "k": 2, "mode": "projection",
                              "project": ["Simulation", "Simulation.1"]},
    "family_supervised": {"prefix": SMALL, "pops": ["French", "Yoruba", "Simulation", "Simulation.1"], "k": 2,
                              "mode": "supervised", "references": ["French", "Yoruba"]},
    "aadr_1kg_k3": {"prefix": AADR, "pops": KGP, "k": 3, "mode": "unsupervised", "chromosomes": "20-22"},
    "aadr_1kg_projection": {"prefix": AADR, "pops": KGP, "k": 3, "mode": "projection", "project": ["GIH.DG"],
                                "chromosomes": "20-22"},
}


def to_bed(prefix: str, pops: list[str], work: Path, name: str, chromosomes: str | None,
           extract: Path | None = None) -> Path:
    (work / f"{name}.pops").write_text("\n".join(pops) + "\n")
    out = work / (name + "_all" if chromosomes or extract else name)
    (work / f"{name}.par").write_text(
        f"genotypename: {prefix}.geno\nsnpname: {prefix}.snp\nindivname: {prefix}.ind\n"
        f"outputformat: PACKEDPED\ngenotypeoutname: {out}.bed\nsnpoutname: {out}.bim\n"
        f"indivoutname: {out}.fam\npoplistname: {work / (name + '.pops')}\n")
    subprocess.run([CONVERTF, "-p", str(work / f"{name}.par")], capture_output=True, check=True)
    if chromosomes or extract:
        # ADMIXTURE refuses SNPs with no genotype at all
        args = ["--chr", chromosomes, "--geno", "0.9999"] if extract is None else ["--extract", str(extract)]
        subprocess.run([PLINK2, "--bfile", str(out), *args, "--make-bed", "--out", str(work / name)],
                       capture_output=True, check=True)
        out = work / name
    return out


def common_snps(work: Path, beds: list[Path]) -> Path:
    """SNPs present in every one of `beds` (each already without all-missing SNPs)."""
    sets = [[line.split()[1] for line in Path(f"{b}.bim").read_text().splitlines()] for b in beds]
    keep = set(sets[0]).intersection(*sets[1:])
    path = work / "common.snps"
    path.write_text("".join(f"{s}\n" for s in sets[0] if s in keep))
    return path


def run_admixture(bed: Path, k: int, extra: list[str]) -> dict:
    r = subprocess.run([ADMIXTURE, *FLAGS, *extra, f"{bed}.bed", str(k)], cwd=bed.parent, capture_output=True,
                       text=True, check=True)
    log = r.stdout
    loglik = float(re.findall(r"^Loglikelihood: (\S+)", log, re.M)[-1])
    iters = re.findall(r"Converged in (\d+) iterations", log)
    cv = re.findall(r"CV error \(K=\d+\): (\S+)", log)
    base = bed.parent / bed.name
    return {
        "loglik": loglik,
        "iterations": int(iters[-1]) if iters else None,
        "cv_error": float(cv[0]) if cv else None,
        "Q": np.loadtxt(f"{base}.{k}.Q", ndmin=2).tolist(),
        "P": np.loadtxt(f"{base}.{k}.P", ndmin=2).tolist(),
    }


def run(name: str) -> None:
    cfg = CONFIGS[name]
    print(f"== {name}", flush=True)
    k = cfg["k"]
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        bed = to_bed(cfg["prefix"], cfg["pops"], work, "data", cfg.get("chromosomes"))
        ids = [line.split()[1] for line in Path(f"{bed}.fam").read_text().splitlines()]
        pop_of = {line.split()[0]: line.split()[2] for line in Path(cfg["prefix"] + ".ind").read_text().splitlines()}
        out = {"snps": None, "admixture_version": "1.3.0", "flags": FLAGS, "prefix": os.path.relpath(cfg["prefix"], REPO)
               if cfg["prefix"].startswith(str(REPO)) else cfg["prefix"], "pops": cfg["pops"], "k": k,
               "mode": cfg["mode"], "chromosomes": cfg.get("chromosomes")}
        if cfg["mode"] == "unsupervised":
            out.update(ids=ids, **run_admixture(bed, k, ["--cv=5"]))
        elif cfg["mode"] == "supervised":
            (work / "data.pop").write_text(
                "".join(f"{pop_of[i] if pop_of[i] in cfg['references'] else '-'}\n" for i in ids))
            out.update(ids=ids, references=cfg["references"], **run_admixture(bed, k, ["--supervised"]))
        else:
            new = to_bed(cfg["prefix"], cfg["project"], work, "new", cfg.get("chromosomes"))
            if cfg.get("chromosomes"):
                # both on the SNPs called in each set of samples
                snps = common_snps(work, [bed, new])
                bed = to_bed(cfg["prefix"], cfg["pops"], work, "data2", None, extract=snps)
                new = to_bed(cfg["prefix"], cfg["project"], work, "new2", None, extract=snps)
                ids = [line.split()[1] for line in Path(f"{bed}.fam").read_text().splitlines()]
            assert Path(f"{new}.bim").read_text() == Path(f"{bed}.bim").read_text()
            fit = run_admixture(bed, k, [])
            shutil.copy(f"{bed}.{k}.P", f"{new}.{k}.P.in")
            proj = run_admixture(new, k, ["-P"])
            out["snps"] = [line.split()[1] for line in Path(f"{new}.bim").read_text().splitlines()]
            out.update(ids=[line.split()[1] for line in Path(f"{new}.fam").read_text().splitlines()],
                       project=cfg["project"], P_in=fit["P"], fit_ids=ids, fit_Q=fit["Q"], fit_loglik=fit["loglik"],
                       loglik=proj["loglik"], iterations=proj["iterations"], Q=proj["Q"])
    dest = HERE / "v1.3.0" / name
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "admixture.json").write_text(json.dumps(out) + "\n")
    print(f"  {len(out['ids'])} samples, loglik {out['loglik']}, cv {out.get('cv_error')}", flush=True)


if __name__ == "__main__":
    for name in sys.argv[1:] or list(CONFIGS):
        run(name)
