"""Golden master test: merge against AdmixTools 8 mergeit.

formats/merge/{a,b}.* hold a SNP on the other strand (with and without the
alleles swapped), an A/T SNP, a third allele, a monomorphic SNP on either
side and two SNPs of file B at one position; formats/merge/reference.* is
mergeit's merge of them (formats/generate_merge.sh). popmix used to drop the
other-strand SNPs and the monomorphic ones without saying why.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import popmix as pm

HERE = Path(__file__).parent / "formats" / "merge"


def _eigenstrat(prefix: Path) -> list[tuple[str, str, str, str]]:
    snps = [line.split() for line in Path(f"{prefix}.snp").read_text().splitlines()]
    geno = Path(f"{prefix}.geno").read_text().splitlines()
    return [(s[0], s[4], s[5], g) for s, g in zip(snps, geno)]


def test_merge_writes_what_mergeit_writes(tmp_path):
    out = tmp_path / "merged"
    result = pm.merge(file_a=str(HERE / "a"), file_b=str(HERE / "b"), output=str(out),
                      output_format="eigenstrat")
    assert _eigenstrat(out) == _eigenstrat(HERE / "reference")
    assert [line.split()[0] for line in Path(f"{out}.ind").read_text().splitlines()] == ["a1", "a2", "b1", "b2"]
    assert (result["strand_ambiguous_snps"], result["snps_mismatch"], result["other_strand_snps"]) == (1, 1, 2)


def test_cli_says_what_was_left_out(tmp_path):
    run = subprocess.run(["popmix", "merge", "-a", str(HERE / "a"), "-b", str(HERE / "b"), "-o", str(tmp_path / "m"),
                          "-of", "eigenstrat"], check=True, capture_output=True, text=True)
    assert "Left out, A/T or C/G (strand unknown): 1" in run.stdout
    assert "Left out, allele mismatch: 1" in run.stdout
    assert "On the other strand (complemented): 2" in run.stdout


def test_merge_streamed_from_packed_and_plink_writes_the_same(tmp_path):
    # packed .geno and PLINK .bed inputs are merged SNP by SNP, read in place
    pm.convertf(str(HERE / "a"), str(tmp_path / "a"), "packedancestrymap")
    pm.convertf(str(HERE / "b"), str(tmp_path / "b"), "plink")
    text = pm.merge(file_a=str(HERE / "a"), file_b=str(HERE / "b"), output=str(tmp_path / "text"),
                    output_format="eigenstrat")
    for fmt in ("eigenstrat", "packedancestrymap", "ped", "tgeno"):
        out = tmp_path / f"streamed_{fmt}"
        streamed = pm.merge(file_a=str(tmp_path / "a"), file_b=str(tmp_path / "b"), output=str(out),
                            output_format=fmt)
        assert streamed == text
        if fmt != "eigenstrat":
            pm.convertf(str(out), str(tmp_path / f"back_{fmt}"), "eigenstrat")
            out = tmp_path / f"back_{fmt}"
        assert _eigenstrat(out) == _eigenstrat(HERE / "reference")


AMBIGUOUS = Path(__file__).parent / "formats" / "merge_ambiguous"


def test_strandcheck_leaves_out_ambiguous_snps_as_mergeit(tmp_path):
    out = tmp_path / "merged"
    pm.merge(file_a=str(AMBIGUOUS / "a"), file_b=str(AMBIGUOUS / "b"), output=str(out), output_format="eigenstrat")
    assert _eigenstrat(out) == _eigenstrat(AMBIGUOUS / "reference_strandcheck_yes")


def test_without_strandcheck_ambiguous_snps_are_kept_as_mergeit(tmp_path):
    out = tmp_path / "merged"
    result = pm.merge(file_a=str(AMBIGUOUS / "a"), file_b=str(AMBIGUOUS / "b"), output=str(out),
                      output_format="eigenstrat", strandcheck=False)
    assert _eigenstrat(out) == _eigenstrat(AMBIGUOUS / "reference_strandcheck_no")
    assert (result["strand_ambiguous_snps"], result["other_strand_snps"], result["snps_mismatch"]) == (0, 0, 1)


def test_the_cli_takes_no_strandcheck(tmp_path):
    out = tmp_path / "m"
    run = subprocess.run(["popmix", "merge", "-a", str(AMBIGUOUS / "a"), "-b", str(AMBIGUOUS / "b"), "-o", str(out),
                          "-of", "eigenstrat", "--no-strandcheck"], check=True, capture_output=True, text=True)
    assert "Left out, allele mismatch: 1" in run.stdout
    assert _eigenstrat(out) == _eigenstrat(AMBIGUOUS / "reference_strandcheck_no")
