"""popmix writes the packed headers as AdmixTools does, byte for byte:
``GENO %7d %7d %x %x`` (and ``TGENO``), so its .geno is the file convertf
writes. tests/golden/dates/family_chr20_22.geno was written by AdmixTools."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

SMALL = "tests/golden/dates/family_chr20_22"
pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def test_packed_is_admixtools_file(tmp_path):
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", str(tmp_path / "p"), "-of", "packedancestrymap"],
                   check=True, capture_output=True)
    assert (tmp_path / "p.geno").read_bytes() == Path(f"{SMALL}.geno").read_bytes()


def test_tgeno_header_is_admixtools_layout(tmp_path):
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", str(tmp_path / "t"), "-of", "tgeno"],
                   check=True, capture_output=True)
    header = (tmp_path / "t.geno").read_bytes()[:48]
    assert header.rstrip(b"\0") == b"TGENO      54   21162 56afe7a2 c3689fb5"
