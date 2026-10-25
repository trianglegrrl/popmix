"""A snp_range outside the file is a ValueError, never a panic or an abort.

read_plink(prefix, (0, 10**9)) used to allocate for 10**9 SNPs and abort the
interpreter.
"""

import subprocess

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
        reader(prefix, snp_range)


@pytest.mark.parametrize("fmt", ["eigenstrat", "plink", "packedancestrymap"])
def test_in_range(prefixes, fmt):
    reader, prefix = prefixes[fmt]
    n = reader(prefix).n_snps
    assert reader(prefix, (1, n)).n_snps == n - 1
