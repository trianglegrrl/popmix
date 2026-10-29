#!/usr/bin/env python3
"""Generate the PLINK 2 reference for popmix's KING-robust kinship.

Usage: python tests/golden/plink2/generate_king.py [config...]

For each config, EIGENSOFT convertf (poplistname) writes the config's
populations as PLINK .bed/.bim/.fam, and PLINK 2 computes KING-robust
kinship two ways:

  --make-king-table counts cols=id,nsnp,hethet,ibs0,ibs1,kinship
      per pair: SNPs with both genotypes called, the counts of het-het,
      opposite-homozygote (IBS0) and het-hom pairs, and the kinship to the
      6 significant digits PLINK prints
  --make-king triangle bin
      the kinship matrix as doubles (full precision)

Writes tests/golden/plink2/<version>/<config>/king.json with the sample IDs
(in .fam order), the pair table in PLINK's row order, and the matrix's lower
triangle in the same order.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
# AADR v66.p1, each individual labelled <population>.<data type>, as
# tests/golden/aadr/build_aadr_v66p1.py builds it in the reference cache; its
# manifest.json (recorded in each AADR reference) has the release checksums.
AADR_DIR = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_v66p1"
AADR = str(AADR_DIR / "aadr_v66p1")
SMALL = str(REPO / "tests/golden/dates/family_chr20_22")
PLINK2 = os.environ.get("PLINK2", "plink2")
CONVERTF = os.environ.get("CONVERTF", "/usr/bin/convertf")

CONFIGS = {
    "family": {"prefix": SMALL, "pops": ["French", "Yoruba", "Simulation", "Simulation.1"]},
    # 1000 Genomes YRI (with relatives) and CEU, and pseudo-haploid ancient
    # samples with missing data (v62.0's Germany_EN_LBK.AG is split by site in
    # v66.p1: its largest part)
    "aadr_yri_ceu_lbk": {"prefix": AADR, "pops": ["YRI.DG", "CEU.DG", "Germany_HalberstadtSonntagsfeld_EN_LBK.AG"]},
}


def plink2_version() -> str:
    out = subprocess.run([PLINK2, "--version"], capture_output=True, text=True, check=True).stdout
    return out.split()[1]


def to_bed(prefix: str, pops: list[str], work: Path) -> Path:
    (work / "pops.txt").write_text("\n".join(pops) + "\n", encoding="utf-8")
    out = work / "data"
    par = work / "convertf.par"
    par.write_text(
        f"genotypename: {prefix}.geno\nsnpname: {prefix}.snp\nindivname: {prefix}.ind\n"
        f"outputformat: PACKEDPED\ngenotypeoutname: {out}.bed\nsnpoutname: {out}.bim\n"
        f"indivoutname: {out}.fam\npoplistname: {work / 'pops.txt'}\n", encoding="utf-8")
    subprocess.run([CONVERTF, "-p", str(par)], capture_output=True, text=True, check=True)
    return out


def run(name: str) -> None:
    cfg = CONFIGS[name]
    print(f"== {name}", flush=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        bed = to_bed(cfg["prefix"], cfg["pops"], work)
        subprocess.run([PLINK2, "--bfile", str(bed), "--make-king-table", "counts",
                        "cols=id,nsnp,hethet,ibs0,ibs1,kinship", "--make-king", "triangle", "bin",
                        "--threads", "4", "--out", str(work / "king")],
                       capture_output=True, text=True, check=True)
        ids = [line.split("\t")[-1] for line in (work / "king.king.id").read_text(encoding="utf-8").splitlines()[1:]]
        table = (work / "king.kin0").read_text(encoding="utf-8").splitlines()
        header = table[0].lstrip("#").split("\t")
        rows = [dict(zip(header, line.split("\t"))) for line in table[1:]]  # noqa: B905 (strict= needs Python 3.10)
        tri = np.fromfile(work / "king.king.bin", dtype=np.float64)
        fam_ids = [line.split()[1] for line in Path(f"{bed}.fam").read_text(encoding="utf-8").splitlines()]
    pop_of = {line.split()[0]: line.split()[2] for line in Path(cfg["prefix"] + ".ind").read_text(encoding="utf-8").splitlines()}
    fam_pops = [pop_of[i] for i in fam_ids]
    n = len(ids)
    assert len(tri) == n * (n - 1) // 2
    assert ids == fam_ids
    out = {
        "plink2_version": plink2_version(),
        "prefix": os.path.relpath(cfg["prefix"], REPO) if cfg["prefix"].startswith(str(REPO)) else cfg["prefix"],
        "pops": cfg["pops"],
        "ids": ids,
        "sample_pops": fam_pops,
        "pairs": {
            "id1": [r["IID1"] for r in rows],
            "id2": [r["IID2"] for r in rows],
            "nsnp": [int(r["NSNP"]) for r in rows],
            "hethet": [int(r["HETHET"]) for r in rows],
            "ibs0": [int(r["IBS0"]) for r in rows],
            "het1_hom2": [int(r["HET1_HOM2"]) for r in rows],
            "het2_hom1": [int(r["HET2_HOM1"]) for r in rows],
            "kinship_printed": [float(r["KINSHIP"]) for r in rows],
        },
        "kinship_lower_triangle": tri.tolist(),
    }
    if cfg["prefix"] == AADR:
        out["dataset"] = json.loads((AADR_DIR / "manifest.json").read_text(encoding="utf-8"))
    dest = HERE / plink2_version() / name
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "king.json").write_text(json.dumps(out) + "\n", encoding="utf-8")
    print(f"  {n} samples, {len(rows)} pairs, max kinship {tri.max():.4f}", flush=True)


if __name__ == "__main__":
    for name in sys.argv[1:] or list(CONFIGS):
        run(name)
