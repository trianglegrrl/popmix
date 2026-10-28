"""Golden master tests: site-frequency spectra vs admixtools2 and dadi.

admixtools2 (tests/golden/admixtools2/v2.0.10/<config>/sfs.json, from
generate_sfs.R): joint_sfs() on packedancestrymap_to_afs(). popmix's
joint_sfs with compat="admixtools2-2.0.10" must give the same table, count
for count; corrected, it also counts the SNPs admixtools2 drops because
their allele count is not exactly an integer in floating point
(ADMIXTOOLS2-0003), and nothing else changes.

dadi 2.4.4 (tests/golden/dadi/v2.4.4/<config>/sfs.json, from
generate_sfs.py): Spectrum.from_data_dict on allele counts from PLINK 2,
projected, unfolded and folded, and S, pi, Watterson's theta and Tajima's D
of the one-population spectra: rtol 1e-10 (masks exact).
"""

import json
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

from .reference_paths import local

ROOT = Path(__file__).parent
REPO = ROOT.parents[1]
A2 = ["family_3pop"] + [pytest.param(c, marks=pytest.mark.slow) for c in ["dg_3pop", "ancient_2pop"]]
DADI = ["family", pytest.param("aadr_1kg", marks=pytest.mark.slow)]


def _prefix(ref: dict) -> str:
    prefix = local(ref["prefix"]) if os.path.isabs(local(ref["prefix"])) else str(REPO / local(ref["prefix"]))
    if not os.path.exists(prefix + ".geno"):
        pytest.skip(f"genotype data not available: {prefix}")
    return prefix


def _a2(name: str) -> dict:
    return json.loads((ROOT / "admixtools2" / "v2.0.10" / name / "sfs.json").read_text())


@pytest.mark.parametrize("name", A2)
def test_joint_sfs_as_admixtools2(name):
    ref = _a2(name)
    ours = pm.joint_sfs(_prefix(ref), ref["pops"], compat="admixtools2-2.0.10")
    for p in ref["pops"]:
        assert ours[p].tolist() == ref["sfs"][p]
    assert ours["n"].tolist() == ref["sfs"]["n"]
    assert ours.attrs["in_sfs"] == ref["n_in_sfs"]
    assert ours.attrs["n_snps"] == ref["n_snps"]


@pytest.mark.parametrize("name", A2)
def test_joint_sfs_corrected_keeps_the_inexact_snps(name):
    """Corrected, every SNP with data in every population is counted."""
    ref = _a2(name)
    prefix = _prefix(ref)
    compat = pm.joint_sfs(prefix, ref["pops"], compat="admixtools2-2.0.10")
    ours = pm.joint_sfs(prefix, ref["pops"])
    assert ours.attrs["not_integer"] == 0
    assert ours.attrs["in_sfs"] == ours.attrs["n_snps"] - ours.attrs["missing"]
    assert compat.attrs["in_sfs"] + compat.attrs["not_integer"] == ours.attrs["in_sfs"]
    assert (ours["n"] >= compat["n"]).all()


def _dadi(name: str) -> dict:
    return json.loads((ROOT / "dadi" / "v2.4.4" / name / "sfs.json").read_text())


@pytest.mark.parametrize("name", DADI)
@pytest.mark.parametrize("folded", [False, True], ids=["unfolded", "folded"])
def test_projected_sfs_as_dadi(name, folded):
    ref = _dadi(name)
    prefix = _prefix(ref)
    for s in ref["spectra"]:
        ours = pm.sfs(prefix, s["pops"], s["projection"], folded=folded)
        want = np.array(s["folded" if folded else "unfolded"]).reshape(ours.data.shape)
        np.testing.assert_allclose(ours.data, want, rtol=1e-10, atol=1e-9,
                                   err_msg=f"{s['pops']} {s['projection']}")
        if folded:
            assert ours.mask.ravel().tolist() == s["folded_mask"]


@pytest.mark.parametrize("name", DADI)
def test_summary_statistics_as_dadi(name):
    ref = _dadi(name)
    prefix = _prefix(ref)
    for s in ref["spectra"]:
        if len(s["pops"]) != 1:
            continue
        ours = pm.sfs(prefix, s["pops"], s["projection"])
        for stat, value in [("S", ours.S()), ("pi", ours.pi()), ("watterson_theta", ours.watterson_theta()),
                            ("tajima_d", ours.tajima_d())]:
            assert value == pytest.approx(s[stat], rel=1e-10), f"{s['pops']} {s['projection']} {stat}"
