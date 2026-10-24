"""Golden master tests: popmix f-statistics vs admixtools2 run on the same data.

Reference outputs live in tests/golden/admixtools2/v<version>/<config>/ and
are produced by generate_reference.R from admixtools2's own f2_from_geno, f2,
f3 and f4. f2.json records every stage of the f2 pipeline, so a mismatch
points at the stage that causes it:

- snp_count: SNPs left after maxmiss and polymorphism filtering
- block_lengths: polymorphic SNPs per block, block by block
- f2_blocks: the per-block f2 for every pair
- f2_estimates: f2() estimate and standard error for every pair

f3_f4.json holds f3() and f4() for every population combination, from the
same f2 blocks (f3_estimates, f4_estimates).

qpwave.json holds qpwave() models on the same f2 blocks, with its stages:
qpwave_f4 (the f4 matrix), qpwave_covariance (its jackknife covariance), and
qpwave_rank_tests (chi-squared and p-value for each rank). The covariance and
chi-squared are held to the standard-error tolerance, p-values to 1e-4.

Tolerances are those in CLAUDE.md: f-statistics rtol 1e-10, standard errors
rtol 1e-6. SNP counts and block lengths must match exactly.
"""

import json
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

REF_ROOT = Path(__file__).parent / "admixtools2" / "v2.0.10"
CONFIGS = ["dg_6pop", "ancient_5pop"]

F_RTOL = 1e-10
F_ATOL = 1e-15
SE_RTOL = 1e-6
P_RTOL = 1e-4


def load_reference(name: str, file: str = "f2.json") -> dict:
    return json.loads((REF_ROOT / name / file).read_text())


_results = {}


def _result_for(name: str):
    if name not in _results:
        ref = load_reference(name)
        if not os.path.exists(f"{ref['prefix']}.geno"):
            pytest.skip(f"genotype data not available: {ref['prefix']}")
        params = ref["params"]
        f2 = pm.extract_f2(
            ref["prefix"],
            ref["pops"],
            auto_only=params["auto_only"],
            block_size=params["blgsize"],
            maxmiss=params["maxmiss"],
            poly_only=params["poly_only"] == "f2",
        )
        _results[name] = (ref, f2)
    return _results[name]


def _check_snp_count(name: str):
    ref, f2 = _result_for(name)
    assert sum(f2.block_lengths) == ref["n_snps_polymorphic"]


def _check_block_lengths(name: str):
    ref, f2 = _result_for(name)
    ours, theirs = list(f2.block_lengths), ref["block_lengths"]
    assert len(ours) == len(theirs), f"{len(ours)} blocks != admixtools2 {len(theirs)}"
    diff = [k for k, (a, b) in enumerate(zip(ours, theirs)) if a != b]
    assert not diff, f"first differing block {diff[0]}: {ours[diff[0]]} != {theirs[diff[0]]}"


def _check_f2_blocks(name: str):
    ref, f2 = _result_for(name)
    mismatches = []
    for pair, theirs in ref["f2_blocks"].items():
        a, b = pair.split("|")
        ours = np.asarray(f2.f2_per_block(a, b))
        theirs = np.asarray(theirs, dtype=float)
        if ours.shape != theirs.shape:
            mismatches.append(f"{pair}: {ours.shape} blocks != {theirs.shape}")
            continue
        bad = ~np.isclose(ours, theirs, rtol=F_RTOL, atol=F_ATOL)
        if bad.any():
            k = int(np.argmax(bad))
            mismatches.append(
                f"{pair}: {int(bad.sum())} blocks differ, first {k}: {ours[k]!r} != {theirs[k]!r}"
            )
    assert not mismatches, "\n".join(mismatches)


def _check_f2_estimates(name: str):
    ref, f2 = _result_for(name)
    mismatches = []
    for row in ref["f2"]:
        got = f2.f2(row["pop1"], row["pop2"])
        if not np.isclose(got["estimate"], row["est"], rtol=F_RTOL, atol=F_ATOL):
            mismatches.append(f"{row['pop1']}|{row['pop2']}: est {got['estimate']!r} != {row['est']!r}")
        if not np.isclose(got["std_err"], row["se"], rtol=SE_RTOL, atol=0):
            mismatches.append(f"{row['pop1']}|{row['pop2']}: se {got['std_err']!r} != {row['se']!r}")
    assert not mismatches, "\n".join(mismatches)


def _check_stat_rows(name: str, stat: str, n_pops: int):
    _, f2 = _result_for(name)
    rows = load_reference(name, "f3_f4.json")[stat]
    compute = getattr(f2, stat)
    mismatches = []
    for row in rows:
        pops = [row[f"pop{i}"] for i in range(1, n_pops + 1)]
        got = compute(*pops)
        label = ",".join(pops)
        if not np.isclose(got["estimate"], row["est"], rtol=F_RTOL, atol=F_ATOL):
            mismatches.append(f"{stat}({label}): est {got['estimate']!r} != {row['est']!r}")
        if not np.isclose(got["std_err"], row["se"], rtol=SE_RTOL, atol=0):
            mismatches.append(f"{stat}({label}): se {got['std_err']!r} != {row['se']!r}")
    assert not mismatches, f"{len(mismatches)} of {2 * len(rows)} differ, first 10:\n" + "\n".join(
        mismatches[:10]
    )


def _check_f3_estimates(name: str):
    _check_stat_rows(name, "f3", 3)


def _check_f4_estimates(name: str):
    _check_stat_rows(name, "f4", 4)


_qpwave = {}


def _qpwave_models(name: str):
    if name not in _qpwave:
        _, f2 = _result_for(name)
        models = load_reference(name, "qpwave.json")["models"]
        _qpwave[name] = [(m, pm.qpwave(f2, m["left"], m["right"])) for m in models]
    return _qpwave[name]


def _model_label(m: dict) -> str:
    return f"left={','.join(m['left'])} right={','.join(m['right'])}"


def _check_qpwave_f4(name: str):
    for m, result in _qpwave_models(name):
        ours, theirs = np.asarray(result.f4_matrix), np.asarray(m["f4_est"], dtype=float)
        assert ours.shape == theirs.shape, _model_label(m)
        assert np.allclose(ours, theirs, rtol=F_RTOL, atol=F_ATOL), (
            f"{_model_label(m)}: f4 {ours.tolist()} != {theirs.tolist()}"
        )


def _check_qpwave_covariance(name: str):
    for m, result in _qpwave_models(name):
        ours, theirs = np.asarray(result.f4_covariance), np.asarray(m["f4_var"], dtype=float)
        assert ours.shape == theirs.shape, _model_label(m)
        rel = np.max(np.abs(ours - theirs) / np.abs(theirs))
        assert rel <= SE_RTOL, f"{_model_label(m)}: covariance max rel err {rel:.2e}"


def _check_qpwave_rank_tests(name: str):
    mismatches = []
    for m, result in _qpwave_models(name):
        ours = {r.rank: r for r in result.rank_tests}
        for row in m["rankdrop"]:
            r = ours[row["f4rank"]]
            label = f"{_model_label(m)} rank {row['f4rank']}"
            if r.df != row["dof"]:
                mismatches.append(f"{label}: df {r.df} != {row['dof']}")
            if not np.isclose(r.chi_sq, row["chisq"], rtol=SE_RTOL, atol=0):
                mismatches.append(f"{label}: chisq {r.chi_sq!r} != {row['chisq']!r}")
            if not np.isclose(r.p_value, row["p"], rtol=P_RTOL, atol=0):
                mismatches.append(f"{label}: p {r.p_value!r} != {row['p']!r}")
    assert not mismatches, "\n".join(mismatches)


CHECKS = {
    "snp_count": _check_snp_count,
    "block_lengths": _check_block_lengths,
    "f2_blocks": _check_f2_blocks,
    "f2_estimates": _check_f2_estimates,
    "f3_estimates": _check_f3_estimates,
    "f4_estimates": _check_f4_estimates,
    "qpwave_f4": _check_qpwave_f4,
    "qpwave_covariance": _check_qpwave_covariance,
    "qpwave_rank_tests": _check_qpwave_rank_tests,
}


# Known gaps, each with the work item that closes it. strict=True makes the
# test fail as soon as the gap is closed, so the mark gets removed.
KNOWN_GAPS = {
    "qpwave_rank_tests": "item 4: qpWave rank fit (opt_A/opt_B and fudge) not yet matched; "
    "the f4 matrix and covariance already match",
}


def _params():
    for name in CONFIGS:
        for check in CHECKS:
            marks = []
            if check in KNOWN_GAPS:
                marks.append(pytest.mark.xfail(reason=KNOWN_GAPS[check], strict=True))
            yield pytest.param(name, check, marks=marks, id=f"{name}-{check}")


@pytest.mark.slow
@pytest.mark.parametrize("name,check", list(_params()))
def test_fstats_admixtools2(name: str, check: str):
    CHECKS[check](name)
