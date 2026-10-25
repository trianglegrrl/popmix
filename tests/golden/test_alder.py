"""Golden master tests: popmix ALDER vs ALDER 1.03.

The reference (tests/golden/alder/v1.03/<config>/alder.json, from
generate_alder.py) is ALDER 1.03 rebuilt to print full precision
(build_alder.sh; checked against the stock build to its printed precision).
Each config holds the main run (two references) and a one-reference run with
each reference; popmix runs with compat="alder-1.03".

Stages:

- ld_corr: the distance to which test LD is correlated with each
  reference's, and the fit start (exact: these are bin edges)
- curve, jackknife_curves: weighted LD per distance bin (rtol 1e-9, and
  atol 1e-11 of the curve's largest value: ALDER computes it by FFT, as
  popmix, whose rounding is relative to the whole curve) and the affine term; bin counts exact, or 1
  above ALDER's, which truncates its FFT count when printing
- fits: decay, amplitudes and their jackknife SEs at each fit start (rtol 1e-6)
- test: z-scores, p-value, decay differences (rtol 1e-6) and the outcome
- mixture_fraction: the one-reference lower bound (rtol 1e-6)
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import popmix._core as core
import pytest

REF_ROOT = Path(__file__).parent / "alder" / "v1.03"
REPO = Path(__file__).parents[2]
CONFIGS = ["family"] + [pytest.param(c, marks=pytest.mark.slow)
                        for c in ["aadr_asw", "aadr_asw_lbk", "aadr_asw_mindis"]]
WLD_RTOL = 1e-9
FIT_RTOL = 1e-6

_results: dict = {}


def _num(x):
    return np.nan if x is None else x


def _run(name: str, key: str, compat: str | None = "alder-1.03"):
    k = (name, key, compat)
    if k not in _results:
        ref = json.loads((REF_ROOT / name / "alder.json").read_text())
        prefix = ref["prefix"] if os.path.isabs(ref["prefix"]) else str(REPO / ref["prefix"])
        if not os.path.exists(prefix + ".geno"):
            pytest.skip(f"genotype data not available: {prefix}")
        r = ref["runs"][key]
        ours = core.alder(prefix, ref["admixpop"], r["refpops"], mindis=ref["params"].get("mindis"), compat=compat)
        _results[k] = (ref, r, ours)
    return _results[k]


def _keys(name: str) -> list:
    path = REF_ROOT / name / "alder.json"
    if not path.exists():
        return ["main"]
    return list(json.loads(path.read_text())["runs"])


def _cases():
    out = []
    for c in CONFIGS:
        name = c.values[0] if hasattr(c, "values") else c
        marks = c.marks if hasattr(c, "marks") else ()
        for key in _keys(name):
            out.append(pytest.param(name, key, marks=marks, id=f"{name}-{key}"))
    return out


CASES = _cases()


@pytest.mark.parametrize(("name", "key"), CASES)
def test_ld_corr(name, key):
    ref, r, ours = _run(name, key)
    if ref["params"].get("mindis") is not None:
        assert all(s is None for s in ours["ld_corr_stops"])
    else:
        want = [r["ld_corr_stops"][p] for p in r["refpops"]]
        assert [100 * s for s in ours["ld_corr_stops"]] == pytest.approx(want, rel=1e-12)
    if r["fit_start"] is not None:
        assert 100 * ours["run"]["fit_start"] == pytest.approx(r["fit_start"], rel=1e-12)


def _check_curve(ours: dict, theirs: list, label: str):
    t = np.array(theirs, dtype=float)
    d = 100 * np.array(ours["d"][1:])
    w = np.array(ours["wld"][1:])
    n = np.array(ours["count"][1:])
    assert len(d) == len(t), label
    np.testing.assert_allclose(d, t[:, 0], rtol=1e-12, err_msg=label)
    # FFT rounding is relative to the curve's largest values, not to each bin
    scale = np.nanmax(np.abs(t[:-1, 1]))
    np.testing.assert_allclose(w, t[:, 1], rtol=WLD_RTOL, atol=1e-11 * scale, err_msg=label)
    diff = np.round(n[:-1]) - t[:-1, 2]
    assert ((diff == 0) | (diff == 1)).all(), f"{label}: counts differ by {diff[(diff != 0) & (diff != 1)][:5]}"
    # ALDER prints counts with %Ld (a long double format) given a long long:
    # the inter-chromosome pair count, past 2^31, comes out wrapped to 32 bits
    last = int(round(n[-1])) - int(t[-1, 2])
    assert last % 2**32 == 0, f"{label}: inter-chromosome pairs {n[-1]} vs {t[-1, 2]}"


@pytest.mark.parametrize(("name", "key"), CASES)
def test_curve(name, key):
    _, r, ours = _run(name, key)
    _check_curve(ours["run"]["curve"], r["curve"], "full")


@pytest.mark.parametrize(("name", "key"), CASES)
def test_jackknife_curves(name, key):
    _, r, ours = _run(name, key)
    labels = ours["chromosomes"]
    assert sorted(labels) == sorted(r["jackknife_curves"])
    for label, curve in zip(labels, ours["run"]["jackknife_curves"]):  # noqa: B905
        _check_curve(curve, r["jackknife_curves"][label], f"without chromosome {label}")


def _check_fits(ours: list, theirs: list, label: str):
    assert len(ours) == len(theirs), label
    for o, t in zip(ours, theirs):  # noqa: B905
        assert 100 * o["mindis"] == pytest.approx(t["mindis_cM"], abs=5e-3), label
        assert o["inter_chrom_affine"] == t["inter_chrom_affine"]
        for var in ("decay", "amp_tot", "amp_exp", "amp_aff"):
            assert o[var]["value"] == pytest.approx(t[var]["value"], rel=FIT_RTOL), f"{label} {var}"
            assert _num(o[var]["se"]) == pytest.approx(_num(t[var]["se"]), rel=FIT_RTOL, nan_ok=True), \
                f"{label} {var} se"
            reps = t[var]["reps"][:-1]  # ALDER lists the full-data value last
            np.testing.assert_allclose(o[var]["reps"], reps, rtol=FIT_RTOL, err_msg=f"{label} {var} reps")


@pytest.mark.parametrize(("name", "key"), CASES)
def test_fits(name, key):
    _, r, ours = _run(name, key)
    main = f"{len(r['refpops'])}-ref {' '.join(r['refpops'])}"
    _check_fits(ours["run"]["fits"], r["fits"][main], main)
    for one, pop in zip(ours["oneref_runs"], r["refpops"]):  # noqa: B905
        _check_fits(one["fits"], r["fits"][f"1-ref {pop}"], f"1-ref {pop}")


@pytest.mark.parametrize("name", list(CONFIGS))
def test_admixture_test(name):
    ref, r, ours = _run(name, "main")
    t, o = r["test"], ours["test"]
    assert o["succeeds"] == t["succeeds"]
    assert o["consistent"] == t["consistent"]
    assert o["z"] == pytest.approx(t["z"], rel=FIT_RTOL)
    assert o["p"] == pytest.approx(t["p"], rel=FIT_RTOL)
    for (z, pct), want in zip(o["decay_diffs"], t["decay_diffs"]):  # noqa: B905
        assert z == pytest.approx(want["z"], rel=FIT_RTOL, abs=1e-12)
        assert pct == pytest.approx(want["percent"], rel=FIT_RTOL)
    one = [[f["decay"]["value"] / f["decay"]["se"], f["amp_exp"]["value"] / f["amp_exp"]["se"]]
           for f in (run["fits"][run["fit_index"]] for run in ours["oneref_runs"])]
    np.testing.assert_allclose(one, t["oneref_z"], rtol=FIT_RTOL)


@pytest.mark.parametrize("name", list(CONFIGS))
def test_mixture_fraction(name):
    ref, _, _ = _run(name, "main")
    for pop in ref["refpops"]:
        _, r, ours = _run(name, f"oneref_{pop}")
        # ALDER prints it as a percentage
        got = [100 * x for x in ours["mixture_fraction"]]
        assert got == pytest.approx(r["mixture_fraction"], rel=FIT_RTOL)


def test_unknown_compat_is_an_error():
    with pytest.raises(ValueError, match="unknown compat"):
        core.alder(str(REPO / "tests/golden/dates/family_chr20_22"), "Simulation.1", ["French"], compat="alder-2")


def test_corrected_affine_term_matches_where_every_snp_has_a_weight():
    _, _, compat = _run("family", "main")
    _, _, corrected = _run("family", "main", compat=None)
    assert math.isclose(compat["run"]["curve"]["wld"][-1], corrected["run"]["curve"]["wld"][-1], rel_tol=0)


@pytest.mark.slow
def test_corrected_affine_term_leaves_out_unweighted_snps():
    """ALDER-0001: with a reference missing at some SNPs, the corrected affine
    term divides by the pairs it sums over, so it is larger."""
    _, _, compat = _run("aadr_asw_lbk", "main")
    _, _, corrected = _run("aadr_asw_lbk", "main", compat=None)
    a, b = compat["run"]["curve"], corrected["run"]["curve"]
    assert b["count"][-1] < a["count"][-1]
    # same sum, fewer pairs
    assert b["wld"][-1] * b["count"][-1] == pytest.approx(a["wld"][-1] * a["count"][-1], rel=1e-12)
