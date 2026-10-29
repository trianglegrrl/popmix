"""The family data are what build_family.py makes from AADR v66.p1.

tests/golden/dates/build_family.py simulates the admixed samples with a
fixed seed; rerun into a temporary directory, it must write byte-identical
files, with the checksums family_truth.json records for the files in the
repository.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from .reference_paths import cache, missing

DATES = Path(__file__).parent / "dates"


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


@pytest.mark.slow
def test_build_is_reproducible(tmp_path):
    if missing(cache() / "aadr_v66p1" / "aadr_v66p1.geno"):
        pytest.skip("AADR v66.p1 subset not built (tests/golden/aadr/build_aadr_v66p1.py)")
    spec = importlib.util.spec_from_file_location("build_family", DATES / "build_family.py")
    build = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build)
    build.main(tmp_path)
    ours = json.loads((tmp_path / "family_truth.json").read_text(encoding="utf-8"))
    shipped = json.loads((DATES / "family_truth.json").read_text(encoding="utf-8"))
    # the numpy that built it is recorded; the files do not depend on it
    # (numpy 1.26.4 built them, 2.0.2 rebuilds them byte for byte)
    ours["simulation"].pop("numpy")
    shipped["simulation"].pop("numpy")
    assert ours == shipped
    for ext in ("geno", "snp", "ind"):
        assert md5(DATES / f"family_chr20_22.{ext}") == shipped["outputs"][f"family_chr20_22.{ext}"]
