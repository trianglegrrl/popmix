"""Golden master tests: popmix READ-style kinship vs READv2.

The reference (tests/golden/readv2/<commit>/<config>/read.json, from
generate_read.py) is READv2 run with its defaults on the config's samples
written as PLINK files by EIGENSOFT convertf. aadr_gurgy is the Neolithic
Gurgy group and other Neolithic France samples (AADR v66.p1 France_Yonne_N.AG
and France_N.AG: 90 pseudo-haploid samples, a pedigree with 78 first-degree
pairs); family is the 49 diploid French and Yoruba samples of the small data.

popmix runs with compat="readv2-2.01" (READV2-0001 to 0003 in DEVIATIONS.md):

- pairs, relationship and first-degree type: exact
- overlap SNPs: exact
- P0, normalized P0, their SEs, z scores, the share of first-degree windows,
  effective SNPs, kinship and its interval: rtol 1e-10
- the normalization value: rtol 1e-12

and corrected (the default): the same P0 and kinship, different SEs and
interval.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

from .reference_paths import local

REF_ROOT = Path(__file__).parent / "readv2" / "73739be"
REPO = Path(__file__).parents[2]
CONFIGS = ["family", pytest.param("aadr_gurgy", marks=pytest.mark.slow)]
RTOL = 1e-10

_results: dict = {}


def _num(values) -> np.ndarray:
    return np.array([np.nan if v is None or (isinstance(v, float) and math.isnan(v)) else float(v)
                     for v in values])


def _run(name: str, compat: str | None = "readv2-2.01"):
    key = (name, compat)
    if key not in _results:
        ref = json.loads((REF_ROOT / name / "read.json").read_text())
        prefix = local(ref["prefix"]) if os.path.isabs(local(ref["prefix"])) else str(REPO / local(ref["prefix"]))
        if not os.path.exists(prefix + ".geno"):
            pytest.skip(f"genotype data not available: {prefix}")
        ours = pm.read_kinship(prefix, ref["pops"], compat=compat)
        # READv2 lists pairs sorted by "id1,id2"
        ours["pair"] = ours["id1"] + "," + ours["id2"]
        ours = ours.set_index("pair").loc[ref["results"]["PairIndividuals"]]
        _results[key] = (ref, ours)
    return _results[key]


@pytest.mark.parametrize("name", CONFIGS)
def test_norm_value(name):
    ref, ours = _run(name)
    assert ours.attrs["norm_value"] == pytest.approx(ref["norm_value"], rel=1e-12)


@pytest.mark.parametrize("name", CONFIGS)
def test_classification(name):
    ref, ours = _run(name)
    res = ref["results"]
    assert ours["relationship"].tolist() == res["Rel"]
    want = ["" if t is None or (isinstance(t, float) and math.isnan(t)) else t for t in res["1st_Type"]]
    assert ours["first_degree_type"].tolist() == want
    assert ours["overlap_snps"].tolist() == [int(x) for x in res["OverlapNSNPs"]]


COLUMNS = {
    "normalized_p0": "P0_mean",
    "p0": "Nonnormalized_P0",
    "p0_se": "Nonnormalized_P0_serr",
    "first_degree_windows_outside": "Perc_Win_1stdeg_P0",
    "effective_snps": "NSNPsXNorm",
    "kinship": "KinshipCoefficient",
    "z_upper": "Zup",
    "z_lower": "Zdown",
}


@pytest.mark.parametrize("name", CONFIGS)
@pytest.mark.parametrize("col", list(COLUMNS))
def test_values(name, col):
    ref, ours = _run(name)
    np.testing.assert_allclose(ours[col].to_numpy(dtype=float), _num(ref["results"][COLUMNS[col]]), rtol=RTOL,
                               atol=0, err_msg=col)


@pytest.mark.parametrize("name", CONFIGS)
def test_interval(name):
    ref, ours = _run(name)
    lo, hi = zip(*(s.split(",") for s in ref["results"]["KinshipCoefficient_95percCI[min,max]"]))  # noqa: B905 (strict= needs Python 3.10)
    np.testing.assert_allclose(ours["kinship_ci_low"], _num(lo), rtol=RTOL)
    np.testing.assert_allclose(ours["kinship_ci_high"], _num(hi), rtol=RTOL)


@pytest.mark.parametrize("name", CONFIGS)
def test_corrected_changes_only_the_errors(name):
    ref, v201 = _run(name)
    _, ours = _run(name, compat=None)
    # the genome-wide P0 does not use the windows
    np.testing.assert_allclose(ours["p0"], v201["p0"], rtol=1e-14)
    np.testing.assert_allclose(ours["kinship"], v201["kinship"], rtol=1e-14)
    assert (ours["p0_se"] != v201["p0_se"]).any()
    # the interval is around the kinship coefficient
    assert ((ours["kinship_ci_low"] < ours["kinship"]) & (ours["kinship"] < ours["kinship_ci_high"])).all()
