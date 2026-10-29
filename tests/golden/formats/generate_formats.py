#!/usr/bin/env python3
"""Reference for reading every genotype format, written by the reference tools.

The same AADR v66.p1 individuals in each format popmix reads, every file
written by AdmixTools 8 convertf (build_convertf.sh) or PLINK 2, never by
popmix:

- packed: convertf's subset of AADR v66.p1 (TGENO) to the populations below,
  as PACKEDANCESTRYMAP; Russia_Samara_EBA_Yamnaya's samples are
  pseudohaploid
- text: convertf, packed -> EIGENSTRAT
- bed: convertf, packed -> PACKEDPED (PLINK .bed/.bim/.fam)
- ped: convertf, packed -> PED
- vcf: PLINK 2, bed -> bgzipped VCF (autosomes)
- fromvcf_bed, fromvcf_text: the reference conversion of that VCF, PLINK 2
  VCF -> .bed, then convertf -> EIGENSTRAT

and admixtools2 f2_from_geno() of the packed files: per-block f2 with 5 cM
blocks (as the files with a genetic map give) and with 2 Mb blocks (a VCF
has no map). The genotype files go to the reference cache, the f2 to
tests/golden/formats/v<admixtools2 version>/reference.json.

Usage: python tests/golden/formats/generate_formats.py
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

HERE = Path(__file__).parent
CACHE = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "formats"
AADR = Path("/references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB")
POPS = ["Russia_Samara_EBA_Yamnaya", "French", "Han", "Mbuti"]
CONVERTF = os.environ.get("POPMIX_CONVERTF",
                          str(Path.home() / ".local/share/popmix-references/admixtools-8.0.2/bin/convertf"))


def sh(*args: str) -> None:
    subprocess.run(list(args), check=True, stdout=subprocess.DEVNULL)


def convertf(src: Path, out: Path, fmt: str, **extra: str) -> None:
    """convertf from an EIGENSTRAT-style or PLINK prefix to `fmt`."""
    if Path(f"{src}.bed").exists():
        names = {"genotypename": f"{src}.bed", "snpname": f"{src}.bim", "indivname": f"{src}.fam"}
    else:
        names = {"genotypename": f"{src}.geno", "snpname": f"{src}.snp", "indivname": f"{src}.ind"}
    ext = {"PACKEDPED": ("bed", "bim", "fam"), "PED": ("ped", "map", "pedind")}.get(fmt, ("geno", "snp", "ind"))
    par = {**names, "outputformat": fmt, "genooutfilename": f"{out}.{ext[0]}",
           "snpoutfilename": f"{out}.{ext[1]}", "indoutfilename": f"{out}.{ext[2]}", **extra}
    p = out.with_suffix(".par")
    p.write_text("".join(f"{k}: {v}\n" for k, v in par.items()), encoding="utf-8")
    # convertf writes a profile (gmon.out) where it runs
    subprocess.run([CONVERTF, "-p", str(p)], check=True, stdout=subprocess.DEVNULL, cwd=out.parent)


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    pops = CACHE / "pops.txt"
    pops.write_text("\n".join(POPS) + "\n", encoding="utf-8")
    packed = CACHE / "packed"
    convertf(AADR, packed, "PACKEDANCESTRYMAP", poplistname=str(pops))
    convertf(packed, CACHE / "text", "EIGENSTRAT")
    convertf(packed, CACHE / "bed", "PACKEDPED")
    convertf(packed, CACHE / "ped", "PED")
    # autosomes only: PLINK 2 will not read chrX back without sex information
    sh("plink2", "--bfile", str(CACHE / "bed"), "--chr", "1-22", "--export", "vcf", "bgz", "id-paste=iid",
       "--allow-extra-chr", "--out", str(CACHE / "vcf"))
    sh("plink2", "--vcf", str(CACHE / "vcf.vcf.gz"), "--make-bed", "--allow-extra-chr",
       "--out", str(CACHE / "fromvcf_bed"))
    convertf(CACHE / "fromvcf_bed", CACHE / "fromvcf_text", "EIGENSTRAT")

    version = subprocess.run(["/usr/bin/Rscript", "-e", 'cat(as.character(packageVersion("admixtools")))'],
                             check=True, capture_output=True, text=True).stdout.strip()
    out = HERE / f"v{version}" / "reference.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    r = f"""
    suppressPackageStartupMessages({{library(admixtools); library(jsonlite)}})
    pops <- c({", ".join(f'"{p}"' for p in POPS)})
    run <- function(blgsize) {{
      f2 <- f2_from_geno("{packed}", pops = pops, blgsize = blgsize, maxmiss = 0, minmaf = 0, maxmaf = 0.5,
                         auto_only = TRUE, adjust_pseudohaploid = TRUE, poly_only = "f2", apply_corr = TRUE,
                         verbose = FALSE)
      pairs <- combn(pops, 2)
      list(block_lengths = as.numeric(gsub("l", "", dimnames(f2)[[3]])),
           f2 = setNames(lapply(seq_len(ncol(pairs)), function(k) sprintf("%.17g", f2[pairs[1, k], pairs[2, k], ])),
                         apply(pairs, 2, paste, collapse = "|")))
    }}
    ind <- read.table("{packed}.ind", colClasses = "character")
    writeLines(toJSON(list(admixtools_version = "{version}", convertf = "AdmixTools 8.0.2", source = "{AADR}",
                           pops = pops,
                           inds = ind$V1, ind_pops = ind$V3, cm = run(0.05), bp = run(2e6)),
                      auto_unbox = TRUE), "{out}")
    """
    subprocess.run(["/usr/bin/Rscript", "-e", r], check=True)
    print("wrote", out)


if __name__ == "__main__":
    main()
