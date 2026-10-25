"""Golden master test: VCF haploid calls and FILTER.

formats/haploid/reference.* (formats/generate_haploid.sh) is the reference
tools' reading of tests/fixtures/small_vcf_haploid.vcf: PLINK 2 reads a
haploid call (male chrX, say) as homozygous and keeps a record whatever its
FILTER; AdmixTools 8 convertf writes the result. popmix read haploid calls
as missing.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

HERE = Path(__file__).parent
VCF = HERE.parent / "fixtures" / "small_vcf_haploid.vcf"
REF = HERE / "formats" / "haploid" / "reference"


def _eigenstrat(prefix: Path) -> list[tuple[str, str, str, str]]:
    snps = [line.split() for line in Path(f"{prefix}.snp").read_text().splitlines()]
    geno = Path(f"{prefix}.geno").read_text().splitlines()
    return [(s[0], s[4], s[5], g) for s, g in zip(snps, geno)]


def test_convertf_writes_what_plink2_and_convertf_write(tmp_path):
    out = tmp_path / "ours"
    subprocess.run(["popmix", "convertf", "-i", str(VCF), "-if", "vcf", "-o", str(out), "-of", "eigenstrat"],
                   check=True, capture_output=True)
    assert _eigenstrat(out) == _eigenstrat(REF)
