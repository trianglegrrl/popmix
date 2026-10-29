"""Golden master tests: every input format, as the reference tools write it.

formats/v<version>/reference.json (formats/generate_formats.py) holds
admixtools2 f2_from_geno() of AADR v66.p1 individuals, including the
pseudohaploid Russia_Samara_EBA_Yamnaya; the same individuals in every
format popmix reads, each written by AdmixTools 8 convertf or PLINK 2, are in
the reference cache.

- f2 from each format must match admixtools2 block by block at the
  f-statistic tolerance: 5 cM blocks for the formats with a genetic map,
  2 Mb blocks for the VCF.
- popmix convertf must write the genotypes the reference tools write for the
  same conversion (PLINK 2 VCF import, then convertf), and convert the
  packed files to the same genotypes as convertf. Binary output is read back
  with convertf, so a reference tool interprets it.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

from .reference_paths import missing

pytestmark = pytest.mark.slow

REF = Path(__file__).parent / "formats" / "v2.0.10" / "reference.json"
CACHE = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "formats"
CONVERTF = os.environ.get("POPMIX_CONVERTF",
                          str(Path.home() / ".local/share/popmix-references/admixtools-8.0.2/bin/convertf"))

AADR_V66 = Path("/references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB")

# format -> (prefix in the cache, blocks, whether the file carries population labels)
FORMATS = {
    "packed": ("packed", "cm", True),
    "eigenstrat": ("text", "cm", True),
    "plink": ("bed", "cm", False),
    "ped": ("ped", "cm", False),
    "vcf": ("vcf.vcf.gz", "bp", False),
}


@pytest.fixture(scope="module")
def ref() -> dict:
    if missing(REF) or missing(CACHE / "packed.geno"):
        pytest.skip("format reference not generated (tests/golden/formats/generate_formats.py)")
    return json.loads(REF.read_text(encoding="utf-8"))


@pytest.mark.parametrize("fmt", list(FORMATS))
def test_f2_from_each_format_matches_admixtools2(ref: dict, fmt: str) -> None:
    name, blocks, labelled = FORMATS[fmt]
    relabel = {} if labelled else {"inds": ref["inds"], "pops": ref["ind_pops"]}
    ours = pm.extract_f2(str(CACHE / name), ref["pops"], **relabel)
    want = ref[blocks]
    assert list(ours.block_lengths) == [int(x) for x in want["block_lengths"]], fmt
    for pair, values in want["f2"].items():
        a, b = pair.split("|")
        np.testing.assert_allclose(ours.f2_per_block(a, b), np.array([float(x) for x in values]),
                                   rtol=1e-10, atol=1e-15, err_msg=f"{fmt} {pair}")


def _eigenstrat(prefix: Path) -> tuple[list[str], dict[str, tuple[str, str, str]]]:
    """Samples, and per SNP id (counted allele, other allele, genotype digits)."""
    samples = [line.split()[0].split(":")[-1] for line in Path(f"{prefix}.ind").read_text(encoding="utf-8").splitlines()]
    snps = [line.split() for line in Path(f"{prefix}.snp").read_text(encoding="utf-8").splitlines()]
    with open(f"{prefix}.geno", encoding="utf-8") as g:
        rows = {s[0]: (s[4], s[5], line.rstrip("\n")) for s, line in zip(snps, g)}
    return samples, rows


_FLIP = str.maketrans("012", "210")


def _assert_same_genotypes(ours: Path, theirs: Path) -> None:
    our_samples, our_rows = _eigenstrat(ours)
    their_samples, their_rows = _eigenstrat(theirs)
    assert our_samples == their_samples
    assert our_rows.keys() == their_rows.keys()
    wrong = []
    for snp, (a, b, geno) in their_rows.items():
        oa, ob, ogeno = our_rows[snp]
        if (oa, ob) == (b, a):  # the same genotypes counted on the other allele
            ogeno = ogeno.translate(_FLIP)
        elif (oa, ob) != (a, b):
            wrong.append((snp, "alleles", (oa, ob), (a, b)))
            continue
        if ogeno != geno:
            wrong.append((snp, "genotypes"))
    assert not wrong, f"{len(wrong)} of {len(their_rows)} SNPs differ, first: {wrong[:3]}"


def _popmix_convertf(src: Path, out: Path, fmt: str, input_format: str | None = None) -> None:
    args = ["popmix", "convertf", "-i", str(src), "-o", str(out), "-of", fmt]
    if input_format:
        args += ["-if", input_format]
    subprocess.run(args, check=True, capture_output=True)


def _to_text_with_convertf(prefix: Path, out: Path) -> None:
    if Path(f"{prefix}.bed").exists():
        names = {"genotypename": f"{prefix}.bed", "snpname": f"{prefix}.bim", "indivname": f"{prefix}.fam"}
    else:
        names = {"genotypename": f"{prefix}.geno", "snpname": f"{prefix}.snp", "indivname": f"{prefix}.ind"}
    par = out.with_suffix(".par")
    par.write_text("".join(f"{k}: {v}\n" for k, v in {
        **names, "outputformat": "EIGENSTRAT", "genooutfilename": f"{out}.geno",
        "snpoutfilename": f"{out}.snp", "indoutfilename": f"{out}.ind"}.items()), encoding="utf-8")
    # convertf writes a profile (gmon.out) where it runs
    subprocess.run([CONVERTF, "-p", str(par)], check=True, capture_output=True, cwd=out.parent)


@pytest.mark.parametrize("fmt", ["eigenstrat", "plink"])
def test_vcf_conversion_matches_plink2_and_convertf(ref: dict, tmp_path: Path, fmt: str) -> None:
    out = tmp_path / "ours"
    _popmix_convertf(CACHE / "vcf.vcf.gz", out, fmt, "vcf")
    if fmt == "plink":
        _to_text_with_convertf(out, tmp_path / "ours_text")
        out = tmp_path / "ours_text"
    _assert_same_genotypes(out, CACHE / "fromvcf_text")


@pytest.mark.parametrize("fmt", ["eigenstrat", "plink", "packedancestrymap", "tgeno"])
def test_packed_conversion_matches_convertf(ref: dict, tmp_path: Path, fmt: str) -> None:
    out = tmp_path / "ours"
    _popmix_convertf(CACHE / "packed", out, fmt)
    if fmt != "eigenstrat":
        _to_text_with_convertf(out, tmp_path / "ours_text")
        out = tmp_path / "ours_text"
    _assert_same_genotypes(out, CACHE / "text")


def test_aadr_v66_is_tgeno() -> None:
    if not Path(f"{AADR_V66}.geno").exists():
        pytest.skip("AADR v66.p1 not available")
    assert pm.dataset_info(str(AADR_V66))["format"] == "tgeno"


def test_f2_straight_from_the_aadr_tgeno_matches_admixtools2(ref: dict) -> None:
    """The packed reference files are convertf's subset of this TGENO file:
    reading only these individuals' records must give admixtools2's f2."""
    if not Path(f"{AADR_V66}.geno").exists():
        pytest.skip("AADR v66.p1 not available")
    ours = pm.extract_f2(str(AADR_V66), ref["pops"])
    want = ref["cm"]
    assert list(ours.block_lengths) == [int(x) for x in want["block_lengths"]]
    for pair, values in want["f2"].items():
        a, b = pair.split("|")
        np.testing.assert_allclose(ours.f2_per_block(a, b), np.array([float(x) for x in values]),
                                   rtol=1e-10, atol=1e-15, err_msg=pair)
