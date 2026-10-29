"""A packed .geno whose header disagrees with the .ind or .snp file is an error.

Deleting one .ind line used to relabel every later sample silently (the
record length is at least 48 bytes, so it did not change).
"""

import subprocess
from pathlib import Path

import pytest

import popmix as pm

FIXTURE = "tests/fixtures/small_eigenstrat"


@pytest.fixture
def packed(tmp_path):
    prefix = tmp_path / "packed"
    subprocess.run(["popmix", "convertf", "-i", FIXTURE, "-o", str(prefix), "-of", "packedancestrymap"],
                   check=True, capture_output=True)
    return prefix


def _drop_first_line(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    path.write_text("".join(lines[1:]), encoding="utf-8")


def test_intact_files_read(packed):
    pops = sorted({line.split()[2] for line in Path(f"{FIXTURE}.ind").read_text(encoding="utf-8").splitlines()})
    n_inds = len(Path(f"{FIXTURE}.ind").read_text(encoding="utf-8").splitlines())
    assert len(pm.read_genotypes(str(packed)).samples) == n_inds
    got = pm.extract_f2(str(packed), pops[:2], poly_only=False).f2(*pops[:2])
    want = pm.extract_f2(FIXTURE, pops[:2], poly_only=False).f2(*pops[:2])
    assert got == want


@pytest.mark.parametrize("ext", ["ind", "snp"])
def test_a_lost_line_is_an_error(packed, ext):
    _drop_first_line(Path(f"{packed}.{ext}"))
    pops = sorted({line.split()[2] for line in Path(f"{FIXTURE}.ind").read_text(encoding="utf-8").splitlines()})
    with pytest.raises(Exception, match="header gives"):
        pm.read_genotypes(str(packed))
    with pytest.raises(Exception, match="header gives"):
        pm.extract_f2(str(packed), pops[:2], poly_only=False)


