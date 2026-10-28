"""extract_f2's SNP filters against admixtools2 2.0.10 (TASKS.md H1).

The references (tests/golden/admixtools2/v2.0.10/snp_filters/<case>.json,
generate_snp_filters.R) are admixtools2's f2_from_geno and fst() on AADR
v66.p1 with transitions, transversions, minmaf/maxmaf, keepsnps and maxmiss:
the block lengths exactly, every f2 block, f2 and FST at CLAUDE.md's
tolerances.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

REF = Path(__file__).parent / "admixtools2" / "v2.0.10" / "snp_filters"
AADR = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_v66p1" \
    / "aadr_v66p1"
CASES = ["transversions", "transitions", "maf", "transversions_maf_maxmiss", "keepsnps", "keepsnps_maf"]
F_RTOL, F_ATOL, SE_RTOL = 1e-10, 1e-15, 1e-6

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not Path(f"{AADR}.geno").exists(), reason=f"AADR v66.p1 build not available: {AADR}"),
]

_cache: dict = {}


def run(case: str):
    if case not in _cache:
        ref = json.loads((REF / f"{case}.json").read_text())
        p = dict(ref["params"])
        if "keepsnps" in p:
            p["keepsnps"] = (REF / p["keepsnps"]).read_text().split()
        _cache[case] = ref, pm.extract_f2(str(AADR), ref["pops"], **p)
    return _cache[case]


@pytest.mark.parametrize("case", CASES)
def test_block_lengths(case):
    ref, f2 = run(case)
    assert f2.block_lengths == ref["block_lengths"]


@pytest.mark.parametrize("case", CASES)
def test_f2_blocks(case):
    ref, f2 = run(case)
    for pair, theirs in ref["f2_blocks"].items():
        a, b = pair.split("|")
        np.testing.assert_allclose(f2.f2_per_block(a, b), theirs, rtol=F_RTOL, atol=F_ATOL, err_msg=pair)


@pytest.mark.parametrize("case", CASES)
def test_f2(case):
    ref, f2 = run(case)
    for row in ref["f2"]:
        got = f2.f2(row["pop1"], row["pop2"])
        assert got["est"] == pytest.approx(row["est"], rel=F_RTOL, abs=F_ATOL)
        assert got["se"] == pytest.approx(row["se"], rel=SE_RTOL)


@pytest.mark.parametrize("case", CASES)
def test_fst(case):
    ref, f2 = run(case)
    for row in ref["fst"]:
        got = f2.fst(row["pop1"], row["pop2"])
        assert got["est"] == pytest.approx(row["est"], rel=F_RTOL, abs=F_ATOL)
        assert got["se"] == pytest.approx(row["se"], rel=SE_RTOL)
