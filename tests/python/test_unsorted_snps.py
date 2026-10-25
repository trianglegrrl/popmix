"""SNPs out of order are an error where jackknife blocks are built.

A block is a run of neighbouring SNPs; with a position going back, blocks
used to mix distant SNPs silently.
"""

import os
import subprocess
from pathlib import Path

import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def _convert(src: str, dst: Path, fmt: str) -> str:
    subprocess.run(["popmix", "convertf", "-i", src, "-o", str(dst), "-of", fmt], check=True,
                   capture_output=True)
    return str(dst)


def _swap_first_two(path: Path) -> None:
    lines = path.read_text().splitlines(keepends=True)
    lines[0], lines[1] = lines[1], lines[0]
    path.write_text("".join(lines))


@pytest.fixture(scope="module")
def unsorted(tmp_path_factory):
    out = tmp_path_factory.mktemp("unsorted")
    text = _convert(SMALL, out / "text", "eigenstrat")
    for ext in ("snp", "geno"):
        _swap_first_two(Path(f"{text}.{ext}"))
    packed = _convert(text, out / "packed", "packedancestrymap")
    return {"eigenstrat": text, "packedancestrymap": packed}


@pytest.mark.parametrize("fmt", ["eigenstrat", "packedancestrymap"])
def test_extract_f2(unsorted, fmt):
    with pytest.raises(ValueError, match="not sorted"):
        pm.extract_f2(unsorted[fmt], POPS)


@pytest.mark.parametrize("fmt", ["eigenstrat", "packedancestrymap"])
def test_statistics_from_genotypes(unsorted, fmt):
    with pytest.raises(ValueError, match="not sorted"):
        pm.qpdstat(unsorted[fmt], *POPS)
