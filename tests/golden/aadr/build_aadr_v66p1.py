#!/usr/bin/env python3
"""The AADR v66.p1 dataset every AADR golden reads, labelled by data type.

The AADR goldens (admixtools2, ADMIXTURE, dadi, DATES, EIGENSOFT, PLINK 2,
READv2, ROLLOFF, and the CLI, Nextflow and reader tests on AADR) run on AADR
v66.p1, 1240K public panel:

    /references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB.{geno,snp,ind}

v66.p1 is the patch release that withdrew the 161 genomes of Jacobs et al.
2019; v62.0 and v66.0 still hold them and must not be used. Its population
labels carry no data-type suffix (`French`), which is on the sample ID
instead (`S_French-1.DG`). The goldens name populations by label and data
type, as `French.DG`, `Chimp.REF` or `Russia_Samara_EBA_Yamnaya.AG`, so
this script relabels every individual `<label>.<last suffix of its sample
ID>` (`French` + `S_French-1.DG` -> `French.DG`; `GLN210.AG.SG` -> `.SG`; a
sample ID without a suffix keeps its label), and AdmixTools 8 convertf
writes the populations below, with those labels, as PACKEDANCESTRYMAP,
which every reference tool reads (v66.p1 itself is TGENO, which only
AdmixTools 8 reads, and converting all 23,089 samples takes hours). Sample
IDs, SNPs (all 1,233,013, in order) and genotypes are v66.p1's own; no
statistic a golden computes depends on the samples left out. A new golden
on other populations adds them to POPULATIONS and reruns this script.

The inputs are checked against the release's CHECKSUMS.normalized.md5
first, and `manifest.json` next to the output records those checksums, the
populations and the output's checksums; each AADR reference records it as
`dataset`. The output is in the reference cache, never in the repository:

    ~/.cache/popmix-references/aadr_v66p1/aadr_v66p1.{geno,snp,ind}

(or under $POPMIX_REFERENCE_CACHE).

Usage: python tests/golden/aadr/build_aadr_v66p1.py
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

RELEASE = Path("/references/aadr/v66.p1")
SOURCE = RELEASE / "v66.p1_1240K.aadr.patch.PUB"
CHECKSUMS = RELEASE / "CHECKSUMS.normalized.md5"
CACHE = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_v66p1"
OUT = CACHE / "aadr_v66p1"
CONVERTF = os.environ.get("POPMIX_CONVERTF",
                          str(Path.home() / ".local/share/popmix-references/admixtools-8.0.2/bin/convertf"))

# Every population an AADR golden uses, as relabelled.
POPULATIONS = [
    # present-day genomes (admixtools2 configs, DATES, the CLI and readers)
    "Chimp.REF", "Mbuti.DG", "Yoruba.DG", "Dinka.DG", "French.DG", "Basque.DG", "Sardinian.DG",
    "Druze.DG", "Han.DG", "Papuan.DG", "Karitiana.DG", "Ulchi.DG",
    # 1000 Genomes (ADMIXTURE, dadi, EIGENSOFT, PLINK 2, ROLLOFF, the ploidy rule)
    "CEU.DG", "TSI.DG", "FIN.DG", "GIH.DG", "YRI.DG", "ASW.DG", "CHB.DG",
    # ancient samples: pseudo-haploid (.AG) and diploid shotgun (.SG)
    "Russia_Samara_EBA_Yamnaya.AG", "Russia_EBA_Kalmykia_Yamnaya.SG", "Russia_EBA_Yamnaya.SG",
    "Germany_HalberstadtSonntagsfeld_EN_LBK.AG", "Czechia_EBA_CordedWare.AG", "England_MIA.AG",
    "France_N.AG", "France_Yonne_N.AG",
]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while chunk := f.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def relabel(line: str) -> str:
    sample, sex, label = line.split()
    suffix = sample.rsplit(".", 1)[1] if "." in sample else ""
    return f"{sample} {sex} {label}.{suffix}" if suffix else f"{sample} {sex} {label}"


def main() -> None:
    listed = dict(reversed(line.split()) for line in CHECKSUMS.read_text().splitlines() if line.strip())
    inputs = {}
    for ext in ("geno", "snp", "ind"):
        path = Path(f"{SOURCE}.{ext}")
        digest = md5(path)
        if digest != listed[path.name]:
            raise SystemExit(f"{path}: md5 {digest} != {listed[path.name]} in {CHECKSUMS}")
        inputs[ext] = {"path": str(path), "md5": digest}

    CACHE.mkdir(parents=True, exist_ok=True)
    relabelled = CACHE / "relabelled.ind"
    relabelled.write_text("".join(relabel(line) + "\n" for line in Path(f"{SOURCE}.ind").read_text().splitlines()))
    labels = {line.split()[2] for line in relabelled.read_text().splitlines()}
    missing = [p for p in POPULATIONS if p not in labels]
    if missing:
        raise SystemExit(f"not in v66.p1: {missing}")
    pops = CACHE / "populations.txt"
    pops.write_text("".join(f"{p}\n" for p in POPULATIONS))
    par = CACHE / "convertf.par"
    par.write_text(f"genotypename: {SOURCE}.geno\nsnpname: {SOURCE}.snp\nindivname: {relabelled}\n"
                   f"outputformat: PACKEDANCESTRYMAP\ngenooutfilename: {OUT}.geno\n"
                   f"snpoutfilename: {OUT}.snp\nindoutfilename: {OUT}.ind\npoplistname: {pops}\n")
    # convertf writes a profile (gmon.out) where it runs
    subprocess.run([CONVERTF, "-p", str(par)], check=True, stdout=subprocess.DEVNULL, cwd=CACHE)
    if Path(f"{OUT}.snp").read_bytes() != Path(f"{SOURCE}.snp").read_bytes():
        raise SystemExit("convertf changed the SNP set")

    manifest = {
        "release": "AADR v66.p1, 1240K public (doi:10.7910/DVN/FFIDCW)",
        "checksums": str(CHECKSUMS),
        "inputs": inputs,
        "relabel": "<label>.<last suffix of sample ID>",
        "populations": POPULATIONS,
        "outputs": {ext: {"path": f"{OUT}.{ext}", "md5": md5(Path(f"{OUT}.{ext}"))} for ext in ("geno", "snp", "ind")},
    }
    (CACHE / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    n = len(Path(f"{OUT}.ind").read_text().splitlines())
    print(f"wrote {OUT}.{{geno,snp,ind}}: {n} samples of {len(POPULATIONS)} populations")


if __name__ == "__main__":
    main()
