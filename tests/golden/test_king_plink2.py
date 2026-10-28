"""Golden master tests: popmix KING-robust kinship vs PLINK 2.

The reference (tests/golden/plink2/<version>/<config>/king.json, from
generate_king.py) is PLINK 2 run on the config's populations written as PLINK
files by EIGENSOFT convertf: --make-king-table counts (per-pair counts, exact)
and --make-king bin (the kinship matrix as doubles).

- pairs: the same pairs in the same order
- counts: nsnp, hethet, ibs0, het1_hom2, het2_hom1, exact. PLINK 2's
  HET1_HOM2 counts SNPs het in IID2 and homozygous in IID1, the reverse of its
  name (PLINK2-0001 in DEVIATIONS.md); popmix's het1_hom2 is het in id1.
- kinship: rtol 1e-10 against the binary matrix; minus infinity where a sample
  has no heterozygous call (the pseudo-haploid samples of aadr_yri_ceu_lbk)
"""

import json
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

from .reference_paths import local

REF_ROOT = Path(__file__).parent / "plink2" / "v2.0.0-a.7LM"
REPO = Path(__file__).parents[2]
CONFIGS = ["family", pytest.param("aadr_yri_ceu_lbk", marks=pytest.mark.slow)]

_results: dict = {}


def _run(name: str):
    if name not in _results:
        ref = json.loads((REF_ROOT / name / "king.json").read_text())
        prefix = local(ref["prefix"]) if os.path.isabs(local(ref["prefix"])) else str(REPO / local(ref["prefix"]))
        if not os.path.exists(prefix + ".geno"):
            pytest.skip(f"genotype data not available: {prefix}")
        _results[name] = (ref, pm.king(prefix, ref["pops"]))
    return _results[name]


@pytest.mark.parametrize("name", CONFIGS)
def test_pairs(name):
    ref, ours = _run(name)
    assert ours["id1"].tolist() == ref["pairs"]["id1"]
    assert ours["id2"].tolist() == ref["pairs"]["id2"]


# popmix column -> PLINK 2 column
COUNTS = {"nsnp": "nsnp", "hethet": "hethet", "ibs0": "ibs0", "het1_hom2": "het2_hom1",
          "het2_hom1": "het1_hom2"}


@pytest.mark.parametrize("name", CONFIGS)
@pytest.mark.parametrize("col", list(COUNTS))
def test_counts(name, col):
    ref, ours = _run(name)
    assert ours[col].tolist() == ref["pairs"][COUNTS[col]]


@pytest.mark.parametrize("name", CONFIGS)
def test_kinship(name):
    ref, ours = _run(name)
    theirs = np.array(ref["kinship_lower_triangle"])
    np.testing.assert_allclose(ours["kinship"], theirs, rtol=1e-10, atol=0)
    # PLINK prints 6 significant digits
    np.testing.assert_allclose(ours["kinship"], ref["pairs"]["kinship_printed"], rtol=1e-5, atol=1e-6)
