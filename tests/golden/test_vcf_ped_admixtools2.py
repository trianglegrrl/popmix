"""Golden master tests: f2 from VCF and PED input vs admixtools2.

vcf_ped/<name>.json (generate_vcf_ped.py) holds admixtools2 f2_from_geno()
of the small family data and of dg_6pop populations of AADR v66.p1 in PLINK
format, with 5 cM blocks and with 2 Mb blocks; the same genotypes as PLINK
2 exports them in VCF and PED are in the reference cache. popmix's f2
blocks from the PED (whose map gives 5 cM blocks) and from the VCF (no map:
2 Mb blocks, as admixtools2 makes without one) must match at the
f-statistic tolerance, and the blocks exactly.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

REF_DIR = Path(__file__).parent / "admixtools2" / "v2.0.10" / "vcf_ped"


def _cases():
    return [pytest.param(p.stem, marks=[pytest.mark.slow] if p.stem.startswith("aadr") else [], id=p.stem)
            for p in sorted(REF_DIR.glob("*.json"))]


def _check(ours, ref: dict, which: str):
    r = ref[which]
    assert list(ours.block_lengths) == [int(x) for x in r["block_lengths"]], which
    for pair, blocks in r["f2"].items():
        a, b = pair.split("|")
        want = np.array([float(x) for x in blocks])
        np.testing.assert_allclose(ours.f2_per_block(a, b), want, rtol=1e-10, atol=1e-15, err_msg=f"{which} {pair}")


@pytest.mark.parametrize("name", _cases())
def test_ped(name):
    ref = json.loads((REF_DIR / f"{name}.json").read_text())
    ped = f"{ref['prefix']}_ped"
    if not os.path.exists(f"{ped}.ped"):
        pytest.skip("PED test data not generated")
    _check(pm.extract_f2(ped, ref["pops"], auto_only=ref["auto_only"]), ref, "cm")


@pytest.mark.parametrize("name", _cases())
def test_vcf(name):
    ref = json.loads((REF_DIR / f"{name}.json").read_text())
    vcf = f"{ref['prefix']}.vcf.gz"
    if not os.path.exists(vcf):
        pytest.skip("VCF test data not generated")
    # a VCF has samples but no populations
    ours = pm.extract_f2(vcf, ref["pops"], auto_only=ref["auto_only"], inds=ref["inds"], pops=ref["ind_pops"])
    _check(ours, ref, "bp")
