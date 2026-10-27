"""A snp_range outside the file is a ValueError, never a panic or an abort.

read_plink(prefix, snp_range=(0, 10**9)) used to allocate for 10**9 SNPs and abort the
interpreter.
"""

import subprocess

import numpy as np
import pytest

import popmix as pm

FIXTURE = "tests/fixtures/small_eigenstrat"


@pytest.fixture(scope="module")
def prefixes(tmp_path_factory):
    out = tmp_path_factory.mktemp("range")
    made = {"eigenstrat": (pm.read_eigenstrat, FIXTURE)}
    for fmt, reader in [("plink", pm.read_plink), ("packedancestrymap", pm.read_genotypes)]:
        prefix = str(out / fmt)
        subprocess.run(["popmix", "convertf", "-i", FIXTURE, "-o", prefix, "-of", fmt], check=True,
                       capture_output=True)
        made[fmt] = (reader, prefix)
    return made


@pytest.mark.parametrize("fmt", ["eigenstrat", "plink", "packedancestrymap"])
@pytest.mark.parametrize("snp_range", [(0, 10**9), (10**9, 10**9 + 2), (3, 2)])
def test_out_of_range(prefixes, fmt, snp_range):
    reader, prefix = prefixes[fmt]
    with pytest.raises(ValueError, match="SNP range"):
        reader(prefix, snp_range=snp_range)


@pytest.mark.parametrize("fmt", ["eigenstrat", "plink", "packedancestrymap"])
def test_in_range(prefixes, fmt):
    reader, prefix = prefixes[fmt]
    n = reader(prefix).n_snps
    assert reader(prefix, snp_range=(1, n)).n_snps == n - 1


def _write_eigenstrat(prefix, chroms):
    """Two populations of two samples, one SNP per chromosome given, 1 cM apart."""
    with open(f"{prefix}.ind", "w") as f:
        f.writelines(f"S{i} U {pop}\n" for i, pop in enumerate(["PopA", "PopA", "PopB", "PopB"]))
    with open(f"{prefix}.snp", "w") as f:
        f.writelines(f"rs{i} {c} {i / 100} {i * 1000 + 1} A G\n" for i, c in enumerate(chroms))
    rows = ["0122", "2101", "1020", "0211", "2012", "1102"]
    with open(f"{prefix}.geno", "w") as f:
        f.writelines(f"{rows[i % len(rows)]}\n" for i in range(len(chroms)))


def _extract(prefix, out, *args):
    return subprocess.run(["popmix", "extract-f2", "-i", prefix, "-o", str(out), *args],
                          capture_output=True, text=True)


def test_extract_f2_snp_range_takes_every_population_by_default(tmp_path):
    prefix = str(tmp_path / "d")
    _write_eigenstrat(prefix, ["1"] * 6)
    r = _extract(prefix, tmp_path / "f2", "--snp-range", "0:6")
    assert r.returncode == 0, r.stderr
    assert pm.read_f2(str(tmp_path / "f2")).populations == ["PopA", "PopB"]


@pytest.mark.parametrize("args", [["--maxmiss", "0.5"], ["--pseudohaploid", "admixtools2"]])
def test_extract_f2_snp_range_refuses_options_it_cannot_apply(tmp_path, args):
    prefix = str(tmp_path / "d")
    _write_eigenstrat(prefix, ["1"] * 6)
    r = _extract(prefix, tmp_path / "f2", "--snp-range", "0:6", *args)
    assert r.returncode != 0
    assert args[0] in r.stderr
    assert not (tmp_path / "f2").exists()


def test_extract_f2_snp_range_refuses_non_autosomal_snps_unless_asked(tmp_path):
    prefix = str(tmp_path / "d")
    _write_eigenstrat(prefix, ["1", "1", "1", "1", "23", "23"])
    r = _extract(prefix, tmp_path / "f2", "--snp-range", "0:6")
    assert r.returncode != 0
    assert "--all-chromosomes" in r.stderr
    assert _extract(prefix, tmp_path / "f2", "--snp-range", "0:4").returncode == 0
    assert _extract(prefix, tmp_path / "f2x", "--snp-range", "0:6", "--all-chromosomes").returncode == 0


@pytest.mark.parametrize("given", ["tests/fixtures/small_vcf.vcf", "tests/fixtures/small_vcf"])
def test_a_snp_range_of_a_vcf_is_those_snps(given):
    whole = pm.read_vcf("tests/fixtures/small_vcf.vcf")
    part = pm.read_genotypes(given, snp_range=(1, 3))
    assert part.n_snps == 2
    assert [s.id for s in part.snps] == [s.id for s in whole.snps][1:3]
    assert np.array_equal(part.to_numpy(), whole.to_numpy()[1:3])
    with pytest.raises(pm.InputError, match="SNP range"):
        pm.read_genotypes(given, snp_range=(3, 9))
