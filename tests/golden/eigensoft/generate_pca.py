#!/usr/bin/env python3
"""Generate the smartpca reference for popmix's PCA.

Usage: python tests/golden/eigensoft/generate_pca.py [config...]

Runs the reference smartpca (EIGENSOFT 7.2.1 printing full precision; see
build_smartpca.sh) and the stock smartpca for each config, checks that the two
agree to the stock binary's printed precision, and writes
tests/golden/eigensoft/v7.2.1/<config>/pca.json with:

  - params: the smartpca parameters (beyond the file names)
  - pca_pops, project_pops: populations the PCA is computed on, and
    populations projected onto it; every other sample is labeled Ignore
  - n_samples_used, n_snps_used: samples and SNPs in the final PCA pass
  - n_snps_killed: SNPs dropped as monomorphic in the samples used
  - eigenvalues: every eigenvalue (evaloutname)
  - samples: each sample in the evec file (PCA samples and projected
    samples; outliers removed are left out), its population and its
    coordinates on the first numoutevec eigenvectors
  - tracy_widom: the TW table smartpca prints (eigenvalue, TW statistic,
    p-value, effective number of markers)
  - outliers: the samples removed as outliers, with iteration, eigenvector
    and score

The relabeled .ind is written to a temporary directory next to symlinks to
the genotype and SNP files; the test recreates it from pca_pops and
project_pops.
"""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
HIPREC = Path(os.environ.get(
    "SMARTPCA_HIPREC",
    Path.home() / ".local/share/popmix-references/eigensoft-7.2.1-hiprec/smartpca"))
STOCK = Path(os.environ.get("SMARTPCA_STOCK", "/usr/lib/eigensoft/smartpca"))
# AADR v66.p1, each individual labelled <population>.<data type>, as
# tests/golden/aadr/build_aadr_v66p1.py builds it in the reference cache; its
# manifest.json (recorded in each AADR reference) has the release checksums.
AADR_DIR = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_v66p1"
AADR = str(AADR_DIR / "aadr_v66p1")
SMALL = str(REPO / "tests/golden/dates/family_chr20_22")

KGP = ["CEU.DG", "TSI.DG", "FIN.DG", "YRI.DG", "CHB.DG", "GIH.DG"]
# v62.0's Germany_EN_LBK.AG is split by site in v66.p1
ANCIENT = ["Germany_HalberstadtSonntagsfeld_EN_LBK.AG", "Czechia_EBA_CordedWare.AG", "England_MIA.AG"]

CONFIGS = {
    # every sample in the PCA, smartpca defaults (5 outlier iterations at 6 sigma)
    "family_all": {"prefix": SMALL, "pca_pops": ["French", "Yoruba", "Simulation", "Simulation.1"],
                       "project_pops": [], "params": {"numoutevec": 4}},
    "family_lsqproject": {"prefix": SMALL, "pca_pops": ["French", "Yoruba"],
                              "project_pops": ["Simulation", "Simulation.1"],
                              "params": {"numoutevec": 4, "lsqproject": "YES"}},
    "family_project_no_outliers": {"prefix": SMALL, "pca_pops": ["French", "Yoruba"],
                                       "project_pops": ["Simulation", "Simulation.1"],
                                       "params": {"numoutevec": 4, "numoutlieriter": 0}},
    "aadr_1kg_lsqproject": {"prefix": AADR, "pca_pops": KGP, "project_pops": ANCIENT,
                                "params": {"numoutevec": 10, "lsqproject": "YES"}},
    "aadr_1kg_project": {"prefix": AADR, "pca_pops": KGP, "project_pops": ANCIENT,
                             "params": {"numoutevec": 10, "numoutlieriter": 0}},
}


def relabel(prefix: str, keep: list[str], out_dir: Path) -> Path:
    """A prefix whose .ind labels every sample outside `keep` Ignore."""
    out = out_dir / Path(prefix).name
    for ext in (".geno", ".snp"):
        os.symlink(os.path.realpath(prefix + ext), str(out) + ext)
    keep_set = set(keep)
    lines = []
    for line in Path(prefix + ".ind").read_text(encoding="utf-8").splitlines():
        ind, sex, pop = line.split()
        lines.append(f"{ind} {sex} {pop if pop in keep_set else 'Ignore'}")
    Path(str(out) + ".ind").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def run_smartpca(binary: Path, prefix: Path, cfg: dict, work: Path, tag: str) -> dict:
    poplist = work / "poplist.txt"
    poplist.write_text("\n".join(cfg["pca_pops"]) + "\n", encoding="utf-8")
    par = {
        "genotypename": f"{prefix}.geno",
        "snpname": f"{prefix}.snp",
        "indivname": f"{prefix}.ind",
        "evecoutname": str(work / f"{tag}.evec"),
        "evaloutname": str(work / f"{tag}.eval"),
        "poplistname": str(poplist),
        **{k: str(v) for k, v in cfg["params"].items()},
    }
    par_file = work / f"{tag}.par"
    par_file.write_text("".join(f"{k}: {v}\n" for k, v in par.items()), encoding="utf-8")
    log = subprocess.run([str(binary), "-p", str(par_file)], capture_output=True, text=True, check=True).stdout
    (work / f"{tag}.log").write_text(log, encoding="utf-8")

    evec_lines = (work / f"{tag}.evec").read_text(encoding="utf-8").splitlines()
    samples = []
    for line in evec_lines[1:]:
        f = line.split()
        samples.append({"id": f[0], "pop": f[-1], "coords": [float(x) for x in f[1:-1]]})
    eigenvalues = [float(x) for x in (work / f"{tag}.eval").read_text(encoding="utf-8").split()]

    tw = []
    in_tw = False
    for line in log.splitlines():
        if line.startswith("## Tracy-Widom"):
            in_tw = True
            continue
        if in_tw:
            f = line.split()
            if f and f[0] == "#N":
                continue
            if not f or not f[0].isdigit():
                break
            tw.append({"n": int(f[0]), "eigenvalue": float(f[1]),
                       "twstat": None if f[3] == "NA" else float(f[3]),
                       "p": None if f[4] == "NA" else float(f[4]),
                       "neff": None if f[5] == "NA" else float(f[5])})
    outliers = [
        {"id": m[1], "iter": int(m[2]), "evec": int(m[3]), "sigmage": float(m[4]), "pop": m[5]}
        for m in re.finditer(r"REMOVED outlier (\S+) iter (\d+) evec (\d+) sigmage (\S+) pop: (\S+)", log)
    ]
    used = re.findall(r"number of samples used: (\d+) number of snps used: (\d+)", log)
    killed = re.findall(r"total number of snps killed in pass: (\d+)\s+used: (\d+)", log)
    return {
        "n_samples_used": len([s for s in samples if s["pop"] in cfg["pca_pops"]]),
        "n_samples_read": int(used[0][0]),
        "n_snps_read": int(used[0][1]),
        "n_snps_killed": int(killed[-1][0]),
        "n_snps_used": int(killed[-1][1]),
        "eigenvalues": eigenvalues,
        "samples": samples,
        "tracy_widom": tw,
        "outliers": outliers,
    }


def check_stock(hiprec: dict, stock: dict) -> None:
    """The stock binary prints coordinates to 4 decimals and eigenvalues to 6."""
    assert [s["id"] for s in hiprec["samples"]] == [s["id"] for s in stock["samples"]]
    for a, b in zip(hiprec["samples"], stock["samples"]):  # noqa: B905 (strict= needs Python 3.10)
        for x, y in zip(a["coords"], b["coords"]):  # noqa: B905 (strict= needs Python 3.10)
            assert abs(x - y) <= 0.5e-4 + 1e-12, (a["id"], x, y)
    for x, y in zip(hiprec["eigenvalues"], stock["eigenvalues"]):  # noqa: B905 (strict= needs Python 3.10)
        assert abs(x - y) <= 0.5e-6 + 1e-12, (x, y)
    assert [o["id"] for o in hiprec["outliers"]] == [o["id"] for o in stock["outliers"]]
    assert hiprec["n_snps_used"] == stock["n_snps_used"]


def main(names: list[str]) -> None:
    for name in names or list(CONFIGS):
        cfg = CONFIGS[name]
        print(f"== {name}", flush=True)
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            prefix = relabel(cfg["prefix"], cfg["pca_pops"] + cfg["project_pops"], work)
            hiprec = run_smartpca(HIPREC, prefix, cfg, work, "hiprec")
            stock = run_smartpca(STOCK, prefix, cfg, work, "stock")
            check_stock(hiprec, stock)
        out = {
            "smartpca_version": "7.2.1",
            "prefix": os.path.relpath(cfg["prefix"], REPO) if cfg["prefix"].startswith(str(REPO)) else cfg["prefix"],
            "pca_pops": cfg["pca_pops"],
            "project_pops": cfg["project_pops"],
            "params": cfg["params"],
            **hiprec,
        }
        if cfg["prefix"] == AADR:
            out["dataset"] = json.loads((AADR_DIR / "manifest.json").read_text(encoding="utf-8"))
        dest = HERE / "v7.2.1" / name
        dest.mkdir(parents=True, exist_ok=True)
        (dest / "pca.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
        print(f"  {out['n_samples_used']} PCA samples, {len(out['samples'])} in evec, "
              f"{out['n_snps_used']} SNPs, {len(out['outliers'])} outliers, "
              f"lambda1 {out['eigenvalues'][0]:.6f}", flush=True)
        assert all(math.isfinite(x) for x in out["eigenvalues"])


if __name__ == "__main__":
    main(sys.argv[1:])
