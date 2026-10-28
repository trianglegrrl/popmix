"""Paths, tools, data subsets and the Case type the benchmark cases share."""

from __future__ import annotations

import os
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
    # an unmeasured run first to warm the page cache; off where one run of
    # the reference takes hours and reading its input seconds
    reference_warmup: bool = True


def _eigenstrat_rows(prefix: Path) -> dict[str, tuple[str, str, str]]:
    snps = [line.split() for line in Path(f"{prefix}.snp").read_text().splitlines()]
    with open(f"{prefix}.geno") as g:
        return {s[0]: (s[4], s[5], line.rstrip("\n")) for s, line in zip(snps, g)}


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


CLASSIC = Path(os.environ.get("POPMIX_ADMIXTOOLS7", HOME / ".local/share/popmix-references/admixtools-7.0.2"))


def _py(script: str, outdir: Path, threads: int, *args: str) -> list[str]:
    path = outdir / "run.py"
    path.write_text(script)
    return ["python", str(path), str(threads), *args, str(outdir)]
