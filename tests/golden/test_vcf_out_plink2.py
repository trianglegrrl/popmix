"""VCF output against PLINK 2 (TASKS.md P1).

formats/vcf_out/reference.vcf is PLINK 2's ``--export vcf`` of
formats/vcf_out/source.* after AdmixTools 8 convertf made it PACKEDPED
(formats/generate_vcf_out.sh), so neither step is popmix's. popmix's VCF
has the same records, GT for GT, but for one field: an unknown allele,
which PLINK 2 writes as its code (X), popmix writes as N, as the VCF
specification asks (DEVIATIONS.md).
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

HERE = Path(__file__).parent / "formats" / "vcf_out"


def records(path: Path) -> list[list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return [line.rstrip("\n").split("\t") for line in f if not line.startswith("##")]


@pytest.fixture(params=["eigenstrat", "packedancestrymap"])
def ours(request, tmp_path) -> Path:
    src = HERE / "source"
    if request.param == "packedancestrymap":
        subprocess.run(["popmix", "convertf", "-i", str(src), "-o", str(tmp_path / "p"), "-of", "packedancestrymap"],
                       check=True, capture_output=True)
        src = tmp_path / "p"
    subprocess.run(["popmix", "convertf", "-i", str(src), "-o", str(tmp_path / "out"), "-of", "vcf"], check=True,
                   capture_output=True)
    return tmp_path / "out.vcf"


def test_records_are_plink2s(ours):
    got, want = records(ours), records(HERE / "reference.vcf")
    assert got[0] == want[0]
    assert len(got) == len(want)
    for g, w in zip(got[1:], want[1:]):
        if w[3] == "X":
            assert g[3] == "N"
            g = g[:3] + ["X"] + g[4:]
        # PLINK 2 marks REF as provisional (INFO PR); popmix writes no INFO
        assert g[:7] + g[8:] == w[:7] + w[8:]


def test_reading_it_back_gives_the_genotypes(ours, tmp_path):
    back = pm.read_genotypes(str(ours))
    source = pm.read_genotypes(str(HERE / "source"))
    # rs5's second allele is unknown: REF N, which a VCF reader of SNPs (popmix's, PLINK 2's
    # --snps-only just-acgt) leaves out
    known = [k for k, snp in enumerate(source.snps) if snp.id != "rs5"]
    assert [x.id for x in back.snps] == [source.snps[k].id for k in known]
    np.testing.assert_array_equal(back.to_numpy(), source.to_numpy()[known])


def test_a_gz_output_is_compressed(tmp_path):
    out = tmp_path / "x.vcf.gz"
    subprocess.run(["popmix", "convertf", "-i", str(HERE / "source"), "-o", str(out), "-of", "vcf"], check=True,
                   capture_output=True)
    assert out.read_bytes()[:2] == b"\x1f\x8b"
    assert records(out) == records(tmp_path / "x.vcf.gz")


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not installed")
def test_bcftools_reads_it(ours):
    r = subprocess.run(["bcftools", "view", "-H", str(ours)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert len(r.stdout.splitlines()) == 8
