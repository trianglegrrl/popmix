"""Golden master tests: popmix f-statistics vs admixtools2 run on the same data.

Reference outputs live in tests/golden/admixtools2/v<version>/<config>/ and
are produced by generate_reference.R from admixtools2's own f2_from_geno, f2,
f3 and f4. f2.json records every stage of the f2 pipeline, so a mismatch
points at the stage that causes it:

- snp_count: SNPs in the f2 blocks after maxmiss and polymorphism filtering
- block_lengths: polymorphic SNPs per block, block by block (with maxmiss > 0,
  after f2_from_geno drops blocks where some pair has no SNP with data)
- f2_blocks: the per-block f2 for every pair
- f2_estimates: f2() estimate and standard error for every pair

f3_f4.json holds f3() and f4() for every population combination, from the
same f2 blocks (f3_estimates, f4_estimates).

qpwave.json holds qpwave() models on the same f2 blocks, with its stages:
qpwave_f4 (the f4 matrix), qpwave_covariance (its jackknife covariance), and
qpwave_rank_tests (chi-squared and p-value for each rank). The covariance and
chi-squared are held to the standard-error tolerance, p-values to 1e-4.

qpadm.json holds qpadm() models, unconstrained and constrained (qpadm_weights:
weights and SEs at the CLAUDE.md weight tolerance, rtol 1e-6; qpadm_fit:
chi-squared and p-value of the model's rank; qpadm_popdrop: the fit of every
subset of the sources), and qpf4ratio() (f4ratio: alpha and SE, rtol 1e-6).
qpadm(constrained = TRUE) fails in admixtools2 with three or more sources, so
for such a model the reference holds only the weights and SEs, computed with
the internals qpadm() uses for them.

fst.json holds fst() for every pair of the config's FST populations (all but
Chimp.REF, a single pseudo-haploid sample for which admixtools2 cannot compute
FST), with its own SNP set and blocks: fst_block_lengths (exact), fst_blocks
(per-block FST, rtol 1e-10) and fst_estimates (est rtol 1e-10, se rtol 1e-6).

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
CONFIGS = ["dg_6pop", "ancient_5pop", "dg_11pop", "ancient_maxmiss"]

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
    assert sum(f2.block_lengths) == sum(ref["block_lengths"])


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
            if "p" not in row:
                # admixtools2 reports p = NA with 0 degrees of freedom; popmix reports 1
                if row["dof"] != 0:
                    mismatches.append(f"{label}: admixtools2 has no p-value")
            elif not np.isclose(r.p_value, row["p"], rtol=P_RTOL, atol=0):
                mismatches.append(f"{label}: p {r.p_value!r} != {row['p']!r}")
    assert not mismatches, "\n".join(mismatches)


# Constrained weights on the boundary are exact in popmix; admixtools2's
# quadprog leaves SEs of ~1e-13 there.
SE_ATOL = 1e-10


def _check_qpadm_weights(name: str):
    _, f2 = _result_for(name)
    mismatches = []
    for m in load_reference(name, "qpadm.json")["models"]:
        r = pm.qpadm(f2, m["target"], m["left"], m["right"], constrained=m["constrained"])
        label = f"target={m['target']} left={','.join(m['left'])} constrained={m['constrained']}"
        for w, se, row in zip(r.weights, r.weight_std_errors, m["weights"]):
            if not np.isclose(w, row["weight"], rtol=SE_RTOL, atol=SE_ATOL):
                mismatches.append(f"{label} {row['left']}: weight {w!r} != {row['weight']!r}")
            if not np.isclose(se, row["se"], rtol=SE_RTOL, atol=SE_ATOL):
                mismatches.append(f"{label} {row['left']}: se {se!r} != {row['se']!r}")
    assert not mismatches, "\n".join(mismatches)


def _check_qpadm_fit(name: str):
    _, f2 = _result_for(name)
    mismatches = []
    for m in load_reference(name, "qpadm.json")["models"]:
        if "rankdrop" not in m:
            continue  # admixtools2 qpadm() fails on this model; weights only
        r = pm.qpadm(f2, m["target"], m["left"], m["right"], constrained=m["constrained"])
        label = f"target={m['target']} left={','.join(m['left'])}"
        row = next(x for x in m["rankdrop"] if x["f4rank"] == len(m["left"]) - 1)
        if r.df != row["dof"]:
            mismatches.append(f"{label}: df {r.df} != {row['dof']}")
        if not np.isclose(r.chi_sq, row["chisq"], rtol=SE_RTOL, atol=0):
            mismatches.append(f"{label}: chisq {r.chi_sq!r} != {row['chisq']!r}")
        if "p" in row and not np.isclose(r.p_value, row["p"], rtol=P_RTOL, atol=0):
            mismatches.append(f"{label}: p {r.p_value!r} != {row['p']!r}")
    assert not mismatches, "\n".join(mismatches)


def _check_qpadm_popdrop(name: str):
    """qpadm()$popdrop: every subset of the sources, as admixtools2 orders it."""
    _, f2 = _result_for(name)
    mismatches = []
    for m in load_reference(name, "qpadm.json")["models"]:
        if "popdrop" not in m:
            continue
        r = pm.qpadm(f2, m["target"], m["left"], m["right"], constrained=m["constrained"])
        label = f"target={m['target']} left={','.join(m['left'])} constrained={m['constrained']}"
        ours = r.popdrop
        if [x.pattern for x in ours] != [row["pat"] for row in m["popdrop"]]:
            mismatches.append(f"{label}: patterns {[x.pattern for x in ours]} != "
                              f"{[row['pat'] for row in m['popdrop']]}")
            continue
        scale = max(abs(row["chisq"]) for row in m["popdrop"])
        for x, row in zip(ours, m["popdrop"]):
            where = f"{label} pat {row['pat']}"
            for field, theirs in (("wt", row["wt"]), ("dof", row["dof"]), ("f4rank", row["f4rank"]),
                                  ("feasible", row["feasible"]), ("best", row.get("best", False))):
                if getattr(x, field) != theirs:
                    mismatches.append(f"{where}: {field} {getattr(x, field)!r} != {theirs!r}")
            if not np.isclose(x.chi_sq, row["chisq"], rtol=SE_RTOL, atol=0):
                mismatches.append(f"{where}: chisq {x.chi_sq!r} != {row['chisq']!r}")
            if "p" in row and not np.isclose(x.p_value, row["p"], rtol=P_RTOL, atol=0):
                mismatches.append(f"{where}: p {x.p_value!r} != {row['p']!r}")
            for source, w in zip(m["left"], x.weights):
                theirs = row.get(source)
                if (w is None) != (theirs is None) or (
                    w is not None and not np.isclose(w, theirs, rtol=SE_RTOL, atol=SE_ATOL)
                ):
                    mismatches.append(f"{where}: weight {source} {w!r} != {theirs!r}")
            if row.get("dofdiff") != x.dofdiff:
                mismatches.append(f"{where}: dofdiff {x.dofdiff!r} != {row.get('dofdiff')!r}")
            for field, tol in (("chisqdiff", dict(rtol=SE_RTOL, atol=SE_RTOL * scale)),
                               ("p_nested", dict(rtol=P_RTOL, atol=0))):
                ours_v, theirs_v = getattr(x, field), row.get(field)
                if (ours_v is None) != (theirs_v is None) or (
                    ours_v is not None and not np.isclose(ours_v, theirs_v, **tol)
                ):
                    mismatches.append(f"{where}: {field} {ours_v!r} != {theirs_v!r}")
    assert not mismatches, "\n".join(mismatches)


def _check_f4ratio(name: str):
    _, f2 = _result_for(name)
    mismatches = []
    for row in load_reference(name, "qpadm.json")["f4ratio"]:
        # admixtools2 qpf4ratio(pops = c(A, O, X, B, C)):
        # alpha = f4(A, O; X, B) / f4(A, O; C, B)
        r = pm.qpf4ratio(f2, test=row["pop3"], source_a=row["pop1"], source_b=row["pop4"],
                         outgroup=row["pop2"], reference=row["pop5"])
        label = ",".join(row[f"pop{i}"] for i in range(1, 6))
        if not np.isclose(r.alpha, row["alpha"], rtol=SE_RTOL, atol=0):
            mismatches.append(f"{label}: alpha {r.alpha!r} != {row['alpha']!r}")
        if not np.isclose(r.std_err, row["se"], rtol=SE_RTOL, atol=0):
            mismatches.append(f"{label}: se {r.std_err!r} != {row['se']!r}")
    assert not mismatches, "\n".join(mismatches)


_fst = {}


def _fst_result_for(name: str):
    if name not in _fst:
        ref = load_reference(name)
        fst_ref = load_reference(name, "fst.json")
        if not os.path.exists(f"{ref['prefix']}.geno"):
            pytest.skip(f"genotype data not available: {ref['prefix']}")
        f2 = pm.extract_f2(
            ref["prefix"],
            fst_ref["pops"],
            block_size=ref["params"]["blgsize"],
            maxmiss=ref["params"]["maxmiss"],
        )
        _fst[name] = (fst_ref, f2)
    return _fst[name]


def _check_fst_block_lengths(name: str):
    ref, f2 = _fst_result_for(name)
    assert list(f2.fst_block_lengths) == ref["block_lengths"]


def _check_fst_blocks(name: str):
    ref, f2 = _fst_result_for(name)
    mismatches = []
    for pair, theirs in ref["fst_blocks"].items():
        ours = np.asarray(f2.fst_per_block(*pair.split("|")))
        theirs = np.asarray(theirs, dtype=float)
        if ours.shape != theirs.shape or not np.allclose(ours, theirs, rtol=F_RTOL, atol=F_ATOL):
            mismatches.append(pair)
    assert not mismatches, f"FST blocks differ for {mismatches}"


def _check_fst_estimates(name: str):
    ref, f2 = _fst_result_for(name)
    mismatches = []
    for row in ref["fst"]:
        got = f2.fst(row["pop1"], row["pop2"])
        label = f"{row['pop1']}|{row['pop2']}"
        if not np.isclose(got["estimate"], row["est"], rtol=F_RTOL, atol=F_ATOL):
            mismatches.append(f"{label}: est {got['estimate']!r} != {row['est']!r}")
        if not np.isclose(got["std_err"], row["se"], rtol=SE_RTOL, atol=0):
            mismatches.append(f"{label}: se {got['std_err']!r} != {row['se']!r}")
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
    "qpadm_weights": _check_qpadm_weights,
    "qpadm_fit": _check_qpadm_fit,
    "qpadm_popdrop": _check_qpadm_popdrop,
    "f4ratio": _check_f4ratio,
    "fst_block_lengths": _check_fst_block_lengths,
    "fst_blocks": _check_fst_blocks,
    "fst_estimates": _check_fst_estimates,
}


# Known gaps, each with the work item that closes it. strict=True makes the
# test fail as soon as the gap is closed, so the mark gets removed.
KNOWN_GAPS: dict = {}


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
