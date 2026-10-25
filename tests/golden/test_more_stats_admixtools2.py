"""Golden master tests: lazadm, qpf4diff and qpfstats vs admixtools2.

more_stats.json (generate_more_stats.R), on the dg_6pop f2 blocks and
genotypes: lazadm (constrained and not), qpf4diff for two sets of five
populations, and qpfstats' fitted f2 blocks for five populations. Weights
at rtol 1e-6, f-statistics at 1e-10, standard errors at 1e-6.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

REF_DIR = Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_6pop"
REF_FILE = REF_DIR / "more_stats.json"
pytestmark = pytest.mark.skipif(not REF_FILE.exists(), reason="more_stats.json not generated")


def _ref() -> dict:
    return json.loads(REF_FILE.read_text())


@pytest.fixture(scope="module")
def f2_blocks():
    f2_ref = json.loads((REF_DIR / "f2.json").read_text())
    if not os.path.exists(f"{f2_ref['prefix']}.geno"):
        pytest.skip(f"genotype data not available: {f2_ref['prefix']}")
    p = f2_ref["params"]
    return pm.extract_f2(f2_ref["prefix"], f2_ref["pops"], auto_only=p["auto_only"], block_size=p["blgsize"],
                         maxmiss=p["maxmiss"], poly_only=p["poly_only"] == "f2")


@pytest.mark.parametrize("k", [0, 1], ids=["constrained", "unconstrained"])
def test_lazadm(f2_blocks, k):
    t = _ref()["lazadm"][k]
    ours = pm.lazadm(f2_blocks, t["left"], t["right"], t["target"], constrained=t["constrained"])
    assert list(ours.left) == [r["left"] for r in t["result"]]
    np.testing.assert_allclose(ours.weight, [float(r["weight"]) for r in t["result"]], rtol=1e-6)
    np.testing.assert_allclose(ours.se, [float(r["se"]) for r in t["result"]], rtol=1e-6)


def test_qpf4diff(f2_blocks):
    t = _ref()["qpf4diff"]
    ours = pm.qpf4diff(f2_blocks, t["pops"])
    np.testing.assert_allclose(ours.f4diff, [float(r["f4diff"]) for r in t["result"]], rtol=1e-10)
    np.testing.assert_allclose(ours.se, [float(r["se"]) for r in t["result"]], rtol=1e-6)
    assert list(ours.n) == [int(r["n"]) for r in t["result"]]
    with pytest.raises(ValueError, match="5 populations"):
        pm.qpf4diff(f2_blocks, ["A", "B"])


@pytest.mark.slow
def test_qpfstats():
    t = _ref()["qpfstats"]
    f2_ref = json.loads((REF_DIR / "f2.json").read_text())
    if not os.path.exists(f"{f2_ref['prefix']}.geno"):
        pytest.skip("genotype data not available")
    ours = pm.qpfstats(f2_ref["prefix"], t["pops"])
    assert list(ours.populations) == t["sorted_pops"]
    assert list(ours.block_lengths) == [int(x) for x in t["block_lengths"]]
    for pair, blocks in t["f2"].items():
        a, b = pair.split("|")
        # NA: a block where some statistic has no SNPs
        want = np.array([np.nan if x == "NA" else float(x) for x in blocks])
        np.testing.assert_allclose(ours.f2_per_block(a, b), want, rtol=1e-9, atol=1e-12 * np.nanmax(np.abs(want)))
