#!/usr/bin/env python3
"""Generate the READv2 reference for popmix's READ-style kinship.

Usage: python tests/golden/readv2/generate_read.py [config...]

READv2 (GuntherLab/READv2, the commit in $READV2/COMMIT) is installed at
$READV2 (default ~/.local/share/popmix-references/readv2) with its own
virtualenv: READ2.py, venv/ with pandas, numpy, matplotlib and plinkio.

For each config, the samples are written as PLINK files by EIGENSOFT
convertf (relabeling first if the config merges populations), READv2 is run
with its defaults (genome-wide P0, 5 Mb jackknife windows, median
normalization), and tests/golden/readv2/<commit>/<config>/read.json records:

  - ids: the samples, in .fam order
  - results: Read_Results.tsv, every column, every pair
  - means: meansP0_AncientDNA_normalized_READv2 (the normalization input)
  - norm_value: the median non-normalized P0 READv2 normalized by
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AADR = "/references/aadr/v62.0_1240k_public"
SMALL = str(REPO / "tests/golden/dates/family_chr20_22")
READV2 = Path(os.environ.get("READV2", Path.home() / ".local/share/popmix-references/readv2"))
CONVERTF = os.environ.get("CONVERTF", "/usr/bin/convertf")


def gurgy_pops() -> list[str]:
    """The pseudo-haploid (.AG) Neolithic France groups: the Gurgy pedigree
    (Rivollat et al. 2023) and its unrelated samples."""
    pops = {line.split()[2] for line in Path(AADR + ".ind").read_text().splitlines()}
    return sorted(p for p in pops if p.startswith("France_N") and p.endswith(".AG") and "_lc" not in p
                  and "contam" not in p and not p.startswith("France_N_BA"))


CONFIGS = {
    "family": {"prefix": SMALL, "pops": ["French", "Yoruba"], "merge": None},
    "aadr_gurgy": {"prefix": AADR, "pops": None, "merge": "France_N_Gurgy"},
}


def to_bed(prefix: str, pops: list[str], merge: str | None, work: Path) -> Path:
    src = prefix
    if merge:
        src = str(work / "relabeled")
        for ext in (".geno", ".snp"):
            os.symlink(os.path.realpath(prefix + ext), src + ext)
        keep = set(pops)
        lines = []
        for line in Path(prefix + ".ind").read_text().splitlines():
            ind, sex, pop = line.split()
            lines.append(f"{ind} {sex} {merge if pop in keep else pop}")
        Path(src + ".ind").write_text("\n".join(lines) + "\n")
        pops = [merge]
    (work / "pops.txt").write_text("\n".join(pops) + "\n")
    out = work / "data"
    (work / "convertf.par").write_text(
        f"genotypename: {src}.geno\nsnpname: {src}.snp\nindivname: {src}.ind\n"
        f"outputformat: PACKEDPED\ngenotypeoutname: {out}.bed\nsnpoutname: {out}.bim\n"
        f"indivoutname: {out}.fam\npoplistname: {work / 'pops.txt'}\n")
    subprocess.run([CONVERTF, "-p", str(work / "convertf.par")], capture_output=True, text=True, check=True)
    return out


def run(name: str) -> None:
    cfg = dict(CONFIGS[name])
    if cfg["pops"] is None:
        cfg["pops"] = gurgy_pops()
    commit = (READV2 / "COMMIT").read_text().strip()
    print(f"== {name}", flush=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        bed = to_bed(cfg["prefix"], cfg["pops"], cfg["merge"], work)
        subprocess.run([str(READV2 / "venv/bin/python"), str(READV2 / "READ2.py"), "-i", "data"], cwd=work,
                       capture_output=True, text=True, check=True)
        results = pd.read_csv(work / "Read_Results.tsv", sep="\t", dtype=str)
        means = pd.read_csv(work / "meansP0_AncientDNA_normalized_READv2", sep="\t", dtype=str)
        ids = [line.split()[1] for line in Path(f"{bed}.fam").read_text().splitlines()]
    nonnorm = means["Nonnormalized_P0"].astype(float)
    out = {
        "read_commit": commit,
        "prefix": os.path.relpath(cfg["prefix"], REPO) if cfg["prefix"].startswith(str(REPO)) else cfg["prefix"],
        "pops": cfg["pops"],
        "merge": cfg["merge"],
        "ids": ids,
        "norm_value": float(nonnorm.median()),
        # strings exactly as READv2 wrote them
        "results": results.to_dict(orient="list"),
        "means": means.to_dict(orient="list"),
    }
    dest = HERE / commit[:7] / name
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "read.json").write_text(json.dumps(out) + "\n")
    print(f"  {len(ids)} samples, {len(results)} pairs:", results["Rel"].value_counts().to_dict(), flush=True)


if __name__ == "__main__":
    for name in sys.argv[1:] or list(CONFIGS):
        run(name)
