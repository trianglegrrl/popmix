"""TGENO (AADR v66 and later): read, detected and written.

tests/fixtures/small_tgeno.geno is AdmixTools 8 convertf's TGENO of
small_eigenstrat; the AADR-scale checks are in tests/golden/test_formats.py.
"""

import subprocess
from pathlib import Path

import numpy as np

import popmix as pm

TGENO = "tests/fixtures/small_tgeno"
TEXT = "tests/fixtures/small_eigenstrat"
SMALL_FAMILY = "tests/golden/dates/family_chr20_22"


def test_detected_and_read():
    assert pm.detect_format(TGENO) == "TGENO"
    assert pm.dataset_info(TGENO)["format"] == "tgeno"
    np.testing.assert_array_equal(pm.read_genotypes(TGENO).to_numpy(), pm.read_genotypes(TEXT).to_numpy())
    assert pm.read_genotypes(TGENO, snp_range=(1, 4)).n_snps == 3
    assert dict(pm.list_populations(TGENO)) == dict(pm.list_populations(TEXT))


def test_f2_and_statistics_match_the_text_file():
    a, b = pm.extract_f2(TGENO, ["PopA", "PopB"], poly_only=False), pm.extract_f2(TEXT, ["PopA", "PopB"],
                                                                                  poly_only=False)
    np.testing.assert_array_equal(a.f2_per_block("PopA", "PopB"), b.f2_per_block("PopA", "PopB"))


def test_convertf_writes_what_convertf_wrote(tmp_path):
    out = tmp_path / "t"
    subprocess.run(["popmix", "convertf", "-i", TEXT, "-o", str(out), "-of", "tgeno"], check=True,
                   capture_output=True)
    ours, theirs = Path(f"{out}.geno").read_bytes(), Path(f"{TGENO}.geno").read_bytes()
    assert ours[48:] == theirs[48:]
    assert ours[:48].rstrip(b"\0").split() == theirs[:48].rstrip(b"\0").split()


def test_dataset_info_lists_populations_in_sample_order(tmp_path):
    """In order of first appearance in the samples, whatever the format, every time."""
    geno = pm.read_genotypes(SMALL_FAMILY)
    expected = list(dict.fromkeys(s.population for s in geno.samples))
    for fmt in ("packedancestrymap", "tgeno", "plink", "eigenstrat"):
        out = str(tmp_path / fmt)
        pm.convertf(SMALL_FAMILY, out, fmt)
        assert pm.dataset_info(out)["populations"] == expected, fmt
