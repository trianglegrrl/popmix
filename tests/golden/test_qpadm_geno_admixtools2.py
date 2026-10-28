"""qpAdm and qpWave on genotypes against admixtools2 2.0.10 (TASKS.md H11).

admixtools2's ``qpadm(prefix)`` computes each f4 from the genotypes, on the
SNPs every population of the model has data for (``allsnps = FALSE``) or on
those of its four populations (``allsnps = TRUE``), rather than from f2
blocks; ``pm.qpadm(prefix, ...)`` does the same. References:
tests/golden/admixtools2/v2.0.10/qpadm_geno/models.json
(generate_qpadm_geno.R), three AADR v66.p1 models with ancient sources.
Weights and SEs at rtol 1e-6, chi-squared at 1e-6, p at 1e-4.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

REF = json.loads((Path(__file__).parent / "admixtools2" / "v2.0.10" / "qpadm_geno" / "models.json").read_text())
AADR = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_v66p1" \
    / "aadr_v66p1"
MODELS = REF["models"]
IDS = [f"{m['target']}-{'allsnps' if m['allsnps'] else 'model'}" for m in MODELS]

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not Path(f"{AADR}.geno").exists(), reason=f"AADR v66.p1 build not available: {AADR}"),
]

_cache: dict = {}


def fit(m: dict) -> pm.QpAdmResult:
    key = (m["target"], m["allsnps"])
    if key not in _cache:
        # admixtools2 takes samples without heterozygous calls in the first 1,000 SNPs as pseudo-haploid
        _cache[key] = pm.qpadm(str(AADR), m["left"], m["right"], m["target"], allsnps=m["allsnps"],
                               pseudohaploid="admixtools2")
    return _cache[key]


@pytest.mark.parametrize("m", MODELS, ids=IDS)
def test_weights_and_ses(m):
    r = fit(m)
    np.testing.assert_allclose(r.weights, m["weights"], rtol=1e-6)
    np.testing.assert_allclose(r.weight_se, m["se"], rtol=1e-6)


@pytest.mark.parametrize("m", MODELS, ids=IDS)
def test_rankdrop(m):
    ours = {t.rank: t for t in fit(m).rankdrop}
    for row in m["rankdrop"]:
        t = ours[row["f4rank"]]
        assert t.dof == row["dof"]
        assert t.chisq == pytest.approx(row["chisq"], rel=1e-6)
        assert t.p == pytest.approx(row["p"], rel=1e-4, abs=1e-300)


@pytest.mark.parametrize("m", MODELS, ids=IDS)
def test_popdrop(m):
    ours = {r.pattern: r for r in fit(m).popdrop}
    for row in m["popdrop"]:
        r = ours[row["pat"]]
        assert r.chisq == pytest.approx(row["chisq"], rel=1e-6), row["pat"]
        assert r.feasible == row["feasible"]


@pytest.mark.parametrize("m", [m for m in MODELS if "qpwave" in m], ids=[i for i, m in zip(IDS, MODELS) if "qpwave" in m])
def test_qpwave(m):
    r = pm.qpwave(str(AADR), [m["target"], *m["left"]], m["right"], pseudohaploid="admixtools2")
    ours = {t.rank: t for t in r.rank_tests}
    for row in m["qpwave"]:
        assert ours[row["f4rank"]].chisq == pytest.approx(row["chisq"], rel=1e-6)


def test_the_counts_are_the_models_blocks_and_snps():
    m = MODELS[0]
    r = fit(m)
    assert r.n_blocks == 713
    assert 1_000_000 < r.n_snps < 1_150_639


def test_allsnps_needs_genotypes():
    f2 = pm.read_admixtools2_f2(Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_6pop" / "f2_dir")
    with pytest.raises(pm.InputError, match="allsnps"):
        pm.qpadm(f2, ["Han.DG", "French.DG"], ["Chimp.REF", "Mbuti.DG", "Yoruba.DG"], "Papuan.DG", allsnps=True)
