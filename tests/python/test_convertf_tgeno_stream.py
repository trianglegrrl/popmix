"""convertf from TGENO reads it in place a part at a time (TASKS.md P4);
its output must be the bytes the same conversion writes from packed
input, for every output format, whole and with a subset."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

SMALL = "tests/golden/dates/family_chr20_22"
FORMATS = ["eigenstrat", "packedancestrymap", "plink", "ancestrymap", "vcf", "ped", "tgeno"]
SUBSET = ["--populations", "Yoruba,Simulation.1", "--chromosomes", "21"]

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def convertf(src: str, out: Path, fmt: str, *extra: str) -> dict[str, bytes]:
    r = subprocess.run(["popmix", "convertf", "-i", src, "-o", str(out), "-of", fmt, *extra],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return {p.name.removeprefix(out.name): p.read_bytes() for p in out.parent.glob(f"{out.name}.*")
            if not p.name.endswith("provenance.json")}


@pytest.fixture(scope="module")
def tgeno(tmp_path_factory) -> str:
    prefix = tmp_path_factory.mktemp("tgeno") / "data"
    convertf(SMALL, prefix, "tgeno")
    return str(prefix)


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("extra", [(), tuple(SUBSET)], ids=["whole", "subset"])
def test_from_tgeno_is_from_packed(tgeno, tmp_path, fmt, extra):
    ours = convertf(tgeno, tmp_path / "from_tgeno", fmt, *extra)
    want = convertf(SMALL, tmp_path / "from_packed", fmt, *extra)
    assert sorted(ours) == sorted(want)
    for ext, data in want.items():
        assert ours[ext] == data, ext
