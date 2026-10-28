"""qpAdm and qpWave at the edges of a model, against admixtools2 2.0.10 (TASKS.md K4).

tests/golden/admixtools2/v2.0.10/degenerate_models/models.json
(generate_degenerate_models.R) holds admixtools2's qpadm() and qpwave() on
AADR v66.p1 f2 blocks for models with no degrees of freedom left (its p is
NA, popmix's NaN), a single source, as many sources as right populations,
and a population in both left and right, which admixtools2 refuses
("Duplicated pops") and popmix refuses too. Weights and SEs at rtol 1e-6,
chi-squared at 1e-6 (and 1e-9 absolute, for the dof-0 fits whose
chi-squared is rounding), p at 1e-4.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

REF = json.loads((Path(__file__).parent / "admixtools2" / "v2.0.10" / "degenerate_models" / "models.json")
                 .read_text())
AADR = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_v66p1" \
    / "aadr_v66p1"
QPADM = {k: m for k, m in REF["qpadm"].items() if "error" not in m}
QPADM_REFUSED = {k: m for k, m in REF["qpadm"].items() if "error" in m}
QPWAVE = {k: m for k, m in REF["qpwave"].items() if "error" not in m}
QPWAVE_REFUSED = {k: m for k, m in REF["qpwave"].items() if "error" in m}

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not Path(f"{AADR}.geno").exists(), reason=f"AADR v66.p1 build not available: {AADR}"),
]


@pytest.fixture(scope="module")
def f2():
    # admixtools2 takes samples without heterozygous calls in the first 1,000 SNPs as pseudo-haploid
    return pm.extract_f2(str(AADR), REF["populations"], pseudohaploid="admixtools2")


def _check_rank_tests(ours, want: list[dict]) -> None:
    by_rank = {t.rank: t for t in ours}
    assert sorted(by_rank) == sorted(row["f4rank"] for row in want)
    for row in want:
        t = by_rank[row["f4rank"]]
        assert t.dof == row["dof"]
        assert t.chisq == pytest.approx(row["chisq"], rel=1e-6, abs=1e-9)
        if row["p"] is None:
            assert math.isnan(t.p)
        else:
            assert t.p == pytest.approx(row["p"], rel=1e-4, abs=1e-300)


@pytest.mark.parametrize("name", list(QPADM))
def test_qpadm(f2, name):
    m = QPADM[name]
    r = pm.qpadm(f2, m["left"], m["right"], m["target"])
    np.testing.assert_allclose(r.weights, np.atleast_1d(m["weights"]), rtol=1e-6)
    np.testing.assert_allclose(r.weight_se, np.atleast_1d(m["se"]), rtol=1e-6, atol=1e-12)
    _check_rank_tests(r.rankdrop, m["rankdrop"])


@pytest.mark.parametrize("name", list(QPWAVE))
def test_qpwave(f2, name):
    m = QPWAVE[name]
    # admixtools2 lists rankdrop from the highest rank down, without the full-rank test
    ours = [t for t in pm.qpwave(f2, m["left"], m["right"]).rank_tests if t.dof > 0]
    _check_rank_tests(ours, m["rankdrop"])


@pytest.mark.parametrize("name", list(QPADM_REFUSED))
def test_qpadm_refuses_what_admixtools2_refuses(f2, name):
    m = QPADM_REFUSED[name]
    with pytest.raises(pm.InputError, match="Duplicate population"):
        pm.qpadm(f2, m["left"], m["right"], m["target"])


@pytest.mark.parametrize("name", list(QPWAVE_REFUSED))
def test_qpwave_refuses_what_admixtools2_refuses(f2, name):
    m = QPWAVE_REFUSED[name]
    with pytest.raises(pm.InputError):
        pm.qpwave(f2, m["left"], m["right"])


def test_the_target_among_the_sources_is_refused(f2):
    # admixtools2 drops it from the sources and fits the rest (DEVIATIONS.md)
    m = QPADM["dof_zero"]
    with pytest.raises(pm.InputError, match="Duplicate population: 'French.DG'"):
        pm.qpadm(f2, [*m["left"], "French.DG"], ["Mbuti.DG", "Han.DG", "Papuan.DG", "Karitiana.DG"], "French.DG")
