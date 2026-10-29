"""Golden master tests: popmix rolloff vs ROLLOFF (ADMIXTOOLS 7.0.2).

The reference (tests/golden/rolloff/v7.0.2/<config>/rolloff.json, from
generate_rolloff.py) is rolloff rebuilt to print full precision
(build_rolloff.sh; checked against the stock build), its curves fit with
rexpfit.r's method (DEoptim start, nls) and its weighted jackknife from
dowtjack, as ADMIXTOOLS' expfit.sh.

- bins: weighted Z-score, regression and correlation per distance bin (rtol
  1e-9) and the pair count (exact), for the full data and each chromosome
  left out
- fit and jackknife (aadr_asw only: with 4 admixed samples the family curve
  determines no date, and nls leaves DEoptim's bounds): popmix's least-squares
  fit (C, A, m) is at least as good as nls's and agrees with it to nls's
  convergence (rtol 1e-4), and the jackknife mean and SE within 1e-4
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

import popmix._core as core

from .reference_paths import local, missing

REF_ROOT = Path(__file__).parent / "rolloff" / "v7.0.2"
REPO = Path(__file__).parents[2]
CONFIGS = ["family", pytest.param("aadr_asw", marks=pytest.mark.slow)]

_results: dict = {}


def _run(name: str):
    if name not in _results:
        ref = json.loads((REF_ROOT / name / "rolloff.json").read_text(encoding="utf-8"))
        prefix = local(ref["prefix"]) if os.path.isabs(local(ref["prefix"])) else str(REPO / local(ref["prefix"]))
        if missing(prefix + ".geno"):
            pytest.skip(f"genotype data not available: {prefix}")
        p = ref["params"]
        ours = core.rolloff(prefix, ref["admixpop"], ref["refs"], binsize=p["binsize"], maxdis=p["maxdis"],
                            runmode=p["runmode"], zdipcorr=p["zdipcorrmode"] == "YES", chithresh=p["chithresh"],
                            chromosomes=ref["chromosomes"])
        _results[name] = (ref, ours)
    return _results[name]


def _check(ours: dict, theirs: list, label: str, first: str):
    t = np.array(theirs)
    np.testing.assert_allclose(100 * np.array(ours["d"]), t[:, 0], rtol=1e-12, err_msg=label)
    np.testing.assert_allclose(ours[first], t[:, 1], rtol=1e-9, atol=1e-12, err_msg=f"{label} {first}")
    np.testing.assert_allclose(ours["regression"], t[:, 2], rtol=1e-9, atol=1e-15, err_msg=f"{label} regression")
    np.testing.assert_allclose(ours["corr"], t[:, 3], rtol=1e-9, atol=1e-15, err_msg=f"{label} corr")
    assert list(ours["count"]) == list(t[:, 4]), label


@pytest.mark.parametrize("name", CONFIGS)
def test_bins(name):
    ref, ours = _run(name)
    _check(ours["bins"], ref["bins"], "full", "zscore")


@pytest.mark.parametrize("name", CONFIGS)
def test_jackknife_bins(name):
    ref, ours = _run(name)
    assert ours["chromosomes"] == list(ref["jackknife_bins"])
    for c, b in zip(ours["chromosomes"], ours["jackknife_bins"]):  # noqa: B905
        # ROLLOFF's jackknife files hold the left-out chromosome's pairs in column 2
        _check(b, ref["jackknife_bins"][c], f"without {c}", "zscore")


def _sse(bins: dict, fit, lval: float = 0.5, hval: float = 100.0) -> float:
    d, y = np.array(bins["d"]) * 100, np.array(bins["corr"])
    m = (d >= lval) & (d <= hval) & np.isfinite(y)
    c, a, k = fit
    return float(np.sum((y[m] - (c + a * np.exp(-k * d[m] / 100))) ** 2))


@pytest.mark.slow
def test_fit_and_jackknife():
    """popmix's fit is at least as good as nls's (its squared error no
    larger, in every case), and the two agree to nls's convergence, 1e-4.

    nls stops once its relative offset falls below 1e-5, which bounds the
    step against the residual, not the parameters: along the flat ridge
    where the affine term, amplitude and decay trade off it stops early. On
    AADR v66.p1 aadr_asw its full fit ends with the decay 1.3e-5 and the
    affine term 4.6e-5 from popmix's, at a squared error 5e-11 above it, and
    the jackknife mean, built from that fit, 6.4e-5 from popmix's."""
    ref, ours = _run("aadr_asw")

    def same(f, r, where):
        assert f == pytest.approx(r, rel=1e-4), where

    assert _sse(ours["bins"], ours["fit"]) <= _sse(ours["bins"], ref["fit"]) * (1 + 1e-12)
    same(ours["fit"], ref["fit"], "all")
    for c, f, b in zip(ours["chromosomes"], ours["jackknife_fits"], ours["jackknife_bins"]):  # noqa: B905
        assert _sse(b, f) <= _sse(b, ref["jackknife_fits"][c]) * (1 + 1e-12), c
        same(f, ref["jackknife_fits"][c], c)
    assert ours["snp_weights"] == [ref["snp_weights"][c] for c in ours["chromosomes"]]
    assert ours["jackknife"] == pytest.approx(ref["jackknife"], rel=1e-4)
