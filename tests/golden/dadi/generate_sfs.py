#!/usr/bin/env python3
"""Generate the dadi reference for popmix's site-frequency spectra.

Usage: python tests/golden/dadi/generate_sfs.py [config...]

The per-population allele counts come from PLINK 2 (`--keep <pop>
--autosome --freq counts` on PLINK files written by EIGENSOFT convertf),
independently of popmix; autosomes only, as PLINK counts males as haploid on
X and Y. PLINK's ALT allele is the allele the .geno counts (popmix's counted
allele), so each SNP enters dadi as calls (REF count, ALT count) with REF as
the ancestral allele.

For each spectrum, dadi (`Spectrum.from_data_dict`, hypergeometric
projection) gives the counts, unfolded and folded, and for one population
its summary statistics (`S`, `pi`, `Watterson_theta`, `Tajima_D`).

Writes tests/golden/dadi/<dadi version>/<config>/sfs.json.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import dadi
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AADR = "/references/aadr/v62.0_1240k_public"
SMALL = str(REPO / "tests/golden/dates/family_chr20_22")
PLINK2 = os.environ.get("PLINK2", "plink2")
CONVERTF = os.environ.get("CONVERTF", "/usr/bin/convertf")

# (populations, projection) per spectrum
CONFIGS = {
    "family": dict(prefix=SMALL, spectra=[
        (["Yoruba"], [42]),
        (["Yoruba"], [20]),
        (["French"], [31]),
        (["French", "Yoruba"], [20, 16]),
    ]),
    "aadr_1kg": dict(prefix=AADR, spectra=[
        (["CEU.DG"], [198]),
        (["CEU.DG"], [40]),
        (["YRI.DG"], [41]),
        (["CEU.DG", "YRI.DG"], [20, 30]),
        (["CEU.DG", "YRI.DG", "CHB.DG"], [10, 12, 8]),
    ]),
}


def counts(prefix: str, pops: list[str], work: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """(REF count, ALT count) per SNP for each population, from PLINK 2."""
    (work / "pops.txt").write_text("\n".join(pops) + "\n")
    bed = work / "data"
    (work / "convertf.par").write_text(
        f"genotypename: {prefix}.geno\nsnpname: {prefix}.snp\nindivname: {prefix}.ind\n"
        f"outputformat: PACKEDPED\ngenotypeoutname: {bed}.bed\nsnpoutname: {bed}.bim\n"
        f"indivoutname: {bed}.fam\npoplistname: {work / 'pops.txt'}\n")
    subprocess.run([CONVERTF, "-p", str(work / "convertf.par")], capture_output=True, check=True)
    pop_of = {line.split()[0]: line.split()[2] for line in Path(prefix + ".ind").read_text().splitlines()}
    fam = [line.split() for line in Path(f"{bed}.fam").read_text().splitlines()]
    out = {}
    for pop in pops:
        keep = work / f"{pop}.keep"
        keep.write_text("".join(f"{f[0]} {f[1]}\n" for f in fam if pop_of[f[1]] == pop))
        subprocess.run([PLINK2, "--bfile", str(bed), "--keep", str(keep), "--autosome", "--freq", "counts",
                        "--out", str(work / pop)], capture_output=True, check=True)
        table = np.loadtxt(work / f"{pop}.acount", dtype=np.int64, comments="#", usecols=(5, 6))
        alt, obs = table[:, 0], table[:, 1]
        out[pop] = (obs - alt, alt)
    return out


def spectrum(c: dict, pops: list[str], projection: list[int], polarized: bool) -> dadi.Spectrum:
    n_snps = len(c[pops[0]][0])
    data = {
        k: {"segregating": ["R", "A"], "outgroup_allele": "R",
            "calls": {p: (int(c[p][0][k]), int(c[p][1][k])) for p in pops}}
        for k in range(n_snps)
    }
    return dadi.Spectrum.from_data_dict(data, pops, projection, polarized=polarized)


def run(name: str) -> None:
    cfg = CONFIGS[name]
    print(f"== {name}", flush=True)
    pops = sorted({p for ps, _ in cfg["spectra"] for p in ps})
    with tempfile.TemporaryDirectory() as tmp:
        c = counts(cfg["prefix"], pops, Path(tmp))
    spectra = []
    for ps, proj in cfg["spectra"]:
        unfolded = spectrum(c, ps, proj, polarized=True)
        folded = spectrum(c, ps, proj, polarized=False)
        entry = {
            "pops": ps,
            "projection": proj,
            "unfolded": np.asarray(unfolded.data).ravel().tolist(),
            "folded": np.asarray(folded.data).ravel().tolist(),
            "folded_mask": np.asarray(folded.mask).ravel().tolist(),
        }
        if len(ps) == 1:
            entry.update(S=float(unfolded.S()), pi=float(unfolded.pi()),
                         watterson_theta=float(unfolded.Watterson_theta()), tajima_d=float(unfolded.Tajima_D()))
        spectra.append(entry)
        print(f"  {ps} -> {proj}: {np.asarray(unfolded.data).sum():.1f} SNPs", flush=True)
    out = {
        "dadi_version": dadi.__version__ if hasattr(dadi, "__version__") else "2.4.4",
        "prefix": os.path.relpath(cfg["prefix"], REPO) if cfg["prefix"].startswith(str(REPO)) else cfg["prefix"],
        "spectra": spectra,
    }
    dest = HERE / "v2.4.4" / name
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "sfs.json").write_text(json.dumps(out) + "\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or list(CONFIGS):
        run(name)
