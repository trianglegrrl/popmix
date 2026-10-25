#!/usr/bin/env python3
"""Reference for f2 from VCF and PED input.

For the small family data and dg_6pop populations of AADR v66.p1 (written
by AdmixTools 8 convertf from tests/golden/aadr/build_aadr_v66p1.py's
dataset, whose manifest the reference records): the genotypes in PLINK
format (popmix convertf), exported by PLINK 2 as VCF and as PED, and
admixtools2 f2_from_geno() of the PLINK files (a copy with
numeric chromosome codes, which admixtools2 needs): with 5 cM blocks
(as the PED's genetic map gives) and with 2 Mb blocks (a VCF has no map,
so popmix, like admixtools2, blocks it by base pairs). The genotype files
go to the reference cache, the f2 blocks to
tests/golden/admixtools2/v2.0.10/vcf_ped/<name>.json.

Usage: python tests/golden/admixtools2/generate_vcf_ped.py
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parents[2]
OUT = HERE / "v2.0.10" / "vcf_ped"
REFERENCES = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references"))
CACHE = REFERENCES / "vcf_ped"
AADR_DIR = REFERENCES / "aadr_v66p1"
CONVERTF = os.environ.get("POPMIX_CONVERTF",
                          str(Path.home() / ".local/share/popmix-references/admixtools-8.0.2/bin/convertf"))
DG6 = ["Han.DG", "French.DG", "Yoruba.DG", "Mbuti.DG", "Papuan.DG", "Chimp.REF"]
DATA = {
    "family": (REPO / "tests/golden/dates/family_chr20_22", ["French", "Yoruba", "Simulation", "Simulation.1"],
               False),
    "aadr_dg6pop": (CACHE / "aadr_dg6pop_packed", DG6, True),
}


def aadr_subset(out: Path, pops: list[str]) -> None:
    """`pops` of the AADR v66.p1 dataset, as PACKEDANCESTRYMAP."""
    poplist = out.with_name(f"{out.name}.pops")
    poplist.write_text("\n".join(pops) + "\n")
    src = AADR_DIR / "aadr_v66p1"
    par = out.with_name(f"{out.name}.par")
    par.write_text(f"genotypename: {src}.geno\nsnpname: {src}.snp\nindivname: {src}.ind\n"
                   f"outputformat: PACKEDANCESTRYMAP\ngenooutfilename: {out}.geno\nsnpoutfilename: {out}.snp\n"
                   f"indoutfilename: {out}.ind\npoplistname: {poplist}\n")
    # convertf writes a profile (gmon.out) where it runs
    subprocess.run([CONVERTF, "-p", str(par)], check=True, stdout=subprocess.DEVNULL, cwd=out.parent)


def sh(*args: str):
    subprocess.run(list(args), check=True, stdout=subprocess.DEVNULL)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    if (AADR_DIR / "aadr_v66p1.geno").exists():
        aadr_subset(DATA["aadr_dg6pop"][0], DG6)
    for name, (src, pops, auto_only) in DATA.items():
        if not Path(f"{src}.geno").exists():
            print("skip", name, "- no", src)
            continue
        bed = CACHE / name
        sh("popmix", "convertf", "-i", str(src), "-o", str(bed), "-of", "plink")
        sh("plink2", "--bfile", str(bed), "--export", "vcf", "bgz", "id-paste=iid", "--allow-extra-chr",
           "--out", str(bed))
        # its own prefix: given a prefix, popmix looks for a VCF first
        sh("plink2", "--bfile", str(bed), "--export", "ped", "--allow-extra-chr", "--out", f"{bed}_ped")
        # admixtools2 cannot parse PLINK's X and Y chromosome codes with
        # auto_only: give it a copy with numeric codes
        num = CACHE / f"{name}_num"
        sh("plink2", "--bfile", str(bed), "--make-bed", "--output-chr", "26", "--allow-extra-chr", "--out", str(num))
        r = f"""
        suppressPackageStartupMessages({{library(admixtools); library(jsonlite)}})
        pops <- c({", ".join(f'"{p}"' for p in pops)})
        run <- function(blgsize) {{
          f2 <- f2_from_geno("{num}", pops = pops, blgsize = blgsize, maxmiss = 0, minmaf = 0, maxmaf = 0.5,
                             auto_only = {"TRUE" if auto_only else "FALSE"}, adjust_pseudohaploid = TRUE,
                             poly_only = "f2", apply_corr = TRUE, verbose = FALSE)
          pairs <- combn(pops, 2)
          list(block_lengths = as.numeric(gsub("l", "", dimnames(f2)[[3]])),
               f2 = setNames(lapply(seq_len(ncol(pairs)), function(k) sprintf("%.17g", f2[pairs[1, k], pairs[2, k], ])),
                             apply(pairs, 2, paste, collapse = "|")))
        }}
        fam <- read.table("{bed}.fam", colClasses = "character")
        writeLines(toJSON(list(pops = pops, auto_only = {"TRUE" if auto_only else "FALSE"}, prefix = "{bed}",
                               inds = fam$V2[fam$V1 %in% pops], ind_pops = fam$V1[fam$V1 %in% pops],
                               cm = run(0.05), bp = run(2e6)), auto_unbox = TRUE), "{OUT / (name + '.json')}")
        """
        subprocess.run(["/usr/bin/Rscript", "-e", r], check=True)
        if name.startswith("aadr"):
            # the f2 values are strings: rewriting the file keeps them exactly
            path = OUT / f"{name}.json"
            ref = json.loads(path.read_text())
            ref["dataset"] = json.loads((AADR_DIR / "manifest.json").read_text())
            path.write_text(json.dumps(ref, separators=(",", ":")) + "\n")
        print("wrote", OUT / f"{name}.json")


if __name__ == "__main__":
    main()
