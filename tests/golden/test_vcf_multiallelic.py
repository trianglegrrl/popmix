"""Golden master test: VCF records that are not biallelic SNPs.

formats/multiallelic/reference.* (formats/generate_multiallelic.sh) is the
reference tools' reading of tests/fixtures/small_vcf_multiallelic.vcf:
PLINK 2 keeps the biallelic SNPs (--snps-only just-acgt --max-alleles 2),
AdmixTools 8 convertf writes them. popmix summed the allele indices of a
multiallelic genotype (0/2 became homozygous ALT) and read indels as
one-base SNPs; it now skips such records, says so, and writes what the
reference tools write.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

import popmix as pm

HERE = Path(__file__).parent
VCF = HERE.parent / "fixtures" / "small_vcf_multiallelic.vcf"
REF = HERE / "formats" / "multiallelic" / "reference"


def _eigenstrat(prefix: Path) -> list[tuple[str, str, str, str]]:
    snps = [line.split() for line in Path(f"{prefix}.snp").read_text().splitlines()]
    geno = Path(f"{prefix}.geno").read_text().splitlines()
    return [(s[0], s[4], s[5], g) for s, g in zip(snps, geno)]


def test_convertf_writes_what_plink2_and_convertf_write(tmp_path):
    out = tmp_path / "ours"
    subprocess.run(["popmix", "convertf", "-i", str(VCF), "-if", "vcf", "-o", str(out), "-of", "eigenstrat"],
                   check=True, capture_output=True)
    assert _eigenstrat(out) == _eigenstrat(REF)


def test_reading_says_what_was_skipped():
    with pytest.warns(UserWarning, match=r"skipped 4 records .*\(1 multiallelic, 3 indels or symbolic"):
        geno = pm.read_genotypes(str(VCF))
    assert geno.n_snps == 3
