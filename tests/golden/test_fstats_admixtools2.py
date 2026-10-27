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
subset of the sources, run with compat="admixtools2-2.0.10";
qpadm_popdrop_corrected: popmix's corrected nested path, ADMIXTOOLS2-0001),
and qpf4ratio() (f4ratio: alpha and SE, rtol 1e-6).
qpadm(constrained = TRUE) fails in admixtools2 with three or more sources
(ADMIXTOOLS2-0002), so for such a model the reference is computed with the
internals qpadm() uses: weights, SEs and popdrop as qpadm() computes them,
and rankdrop from drop_ranks() run unconstrained (rankdrop never uses the
weights).

qpgraph.json holds qpgraph() fits of fixed graphs (seed 1): qpgraph_f3 (the
f3 statistics it fits, est and se, and their inverse covariance ppinv);
qpgraph_fixed (edge lengths, fitted f3 and score at fixed admixture weights,
admixtools2's best and all 0.3, at rtol 1e-10: no optimizer is involved);
qpgraph_fit (with weights optimized, from admixtools2's own random starts:
every start's score at rtol 1e-10 and weights at 1e-6).

fst.json holds fst() for every pair of the config's FST populations (all but
Chimp.REF, a single pseudo-haploid sample for which admixtools2 cannot compute
FST), with its own SNP set and blocks: fst_block_lengths (exact), fst_blocks
(per-block FST, rtol 1e-10) and fst_estimates (est rtol 1e-10, se rtol 1e-6).

popmix reads the genotypes with admixtools2's ploidy rule
(pseudohaploid="admixtools2", ADMIXTOOLS2-0007): in AADR v66.p1 dg_11pop's
Druze.DG holds HGDP00581.DG, a diploid genome with no heterozygous call
among its first 1,000 SNPs, which admixtools2 counts as pseudo-haploid. In
the other configs the two rules agree on every sample.

Tolerances are those in CLAUDE.md: f-statistics rtol 1e-10, standard errors
rtol 1e-6. SNP counts and block lengths must match exactly.
"""

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

REF_ROOT = Path(__file__).parent / "admixtools2" / "v2.0.10"
CONFIGS = ["dg_6pop", "ancient_5pop", "dg_11pop", "ancient_maxmiss", "family_4pop"]

F_RTOL = 1e-10
F_ATOL = 1e-15
SE_RTOL = 1e-6
# Admixture weights (CLAUDE.md)
WEIGHT_TOL = 1e-6
P_RTOL = 1e-4
# A qpgraph fit whose line search failed (L-BFGS-B convergence 52) stopped
# where rounding put it, so its starts agree only this far (see _check_qpgraph_fit)
STALL_RTOL = 1e-5
STALL_WEIGHT_TOL = 1e-3


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
            # admixtools2's ploidy rule (ADMIXTOOLS2-0007): dg_11pop's Druze.DG
            # holds a diploid genome it takes as pseudo-haploid
            pseudohaploid="admixtools2",
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
        r = pm.qpadm(f2, m["left"], m["right"], m["target"], constrained=m["constrained"])
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
            continue
        r = pm.qpadm(f2, m["left"], m["right"], m["target"], constrained=m["constrained"])
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
    """qpadm()$popdrop: every subset of the sources, as admixtools2 orders it.

    Run with compat="admixtools2-2.0.10", which reproduces admixtools2's
    nested path (ADMIXTOOLS2-0001).
    """
    _, f2 = _result_for(name)
    mismatches = []
    for m in load_reference(name, "qpadm.json")["models"]:
        if "popdrop" not in m:
            continue
        r = pm.qpadm(f2, m["left"], m["right"], m["target"], constrained=m["constrained"],
                     compat="admixtools2-2.0.10")
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
            for field, tol in (("chisqdiff", {"rtol": SE_RTOL, "atol": SE_RTOL * scale}),
                               ("p_nested", {"rtol": P_RTOL, "atol": 0})):
                ours_v, theirs_v = getattr(x, field), row.get(field)
                if (ours_v is None) != (theirs_v is None) or (
                    ours_v is not None and not np.isclose(ours_v, theirs_v, **tol)
                ):
                    mismatches.append(f"{where}: {field} {ours_v!r} != {theirs_v!r}")
    assert not mismatches, "\n".join(mismatches)


def _reference_nested_path(rows: list) -> list:
    """The corrected nested path, walked over admixtools2's own model fits.

    From the full model, repeatedly drop the one further source whose model
    is feasible and fits best, down to a single source.
    """
    by_pat = {row["pat"]: row for row in rows}
    path = ["0" * len(rows[0]["pat"])]
    while path[-1].count("0") > 1:
        cur = path[-1]
        children = [cur[:i] + "1" + cur[i + 1:] for i, c in enumerate(cur) if c == "0"]
        children = [c for c in children if by_pat[c]["feasible"]]
        if not children:
            break
        path.append(min(children, key=lambda c: by_pat[c]["chisq"]))
    return path


def _check_qpadm_popdrop_corrected(name: str):
    """ADMIXTOOLS2-0001, live: by default only the nested-path columns differ
    from admixtools2 (isolation), and they follow the nested path through
    admixtools2's own model fits (oracle), each step tested against the model
    it came from with 1 degree of freedom."""
    _, f2 = _result_for(name)
    mismatches = []
    for m in load_reference(name, "qpadm.json")["models"]:
        if "popdrop" not in m:
            continue
        r = pm.qpadm(f2, m["left"], m["right"], m["target"], constrained=m["constrained"])
        compat = pm.qpadm(f2, m["left"], m["right"], m["target"], constrained=m["constrained"],
                          compat="admixtools2-2.0.10")
        label = f"target={m['target']} left={','.join(m['left'])} constrained={m['constrained']}"
        for x, y in zip(r.popdrop, compat.popdrop):
            for field in ("pattern", "wt", "dof", "chi_sq", "p_value", "f4rank", "weights",
                          "feasible"):
                a, b = getattr(x, field), getattr(y, field)
                both_na = isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b)
                if a != b and not both_na:  # a saturated model's p is NA in both
                    mismatches.append(f"{label} pat {x.pattern}: {field} differs from compat")
        rows = {row["pat"]: row for row in m["popdrop"]}
        path = _reference_nested_path(m["popdrop"])
        scale = max(abs(row["chisq"]) for row in m["popdrop"])
        for x in r.popdrop:
            where = f"{label} pat {x.pattern}"
            if x.best != (x.pattern in path):
                mismatches.append(f"{where}: best {x.best} but path is {path}")
            k = path.index(x.pattern) if x.pattern in path else 0
            if k == 0:
                if (x.dofdiff, x.chisqdiff, x.p_nested) != (None, None, None):
                    mismatches.append(f"{where}: has a nested comparison, path is {path}")
                continue
            child, parent = rows[x.pattern], rows[path[k - 1]]
            chisqdiff = child["chisq"] - parent["chisq"]
            p_nested = math.erfc(math.sqrt(chisqdiff / 2)) if chisqdiff > 0 else 1.0
            if x.dofdiff != child["dof"] - parent["dof"] or x.dofdiff != 1:
                mismatches.append(f"{where}: dofdiff {x.dofdiff}")
            if not np.isclose(x.chisqdiff, chisqdiff, rtol=SE_RTOL, atol=SE_RTOL * scale):
                mismatches.append(f"{where}: chisqdiff {x.chisqdiff!r} != {chisqdiff!r}")
            if not np.isclose(x.p_nested, p_nested, rtol=P_RTOL, atol=0):
                mismatches.append(f"{where}: p_nested {x.p_nested!r} != {p_nested!r}")
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


def _qpgraph_reference(name: str) -> list:
    path = REF_ROOT / name / "qpgraph.json"
    if not path.exists():
        pytest.skip(f"{name} has no qpgraph graphs")
    return json.loads(path.read_text())["graphs"]


def _as_list(x) -> list:
    return x if isinstance(x, list) else [x]


def _check_qpgraph_f3(name: str):
    """qpgraph()$f3 est and se, and the inverse f3 covariance it fits with."""
    _, f2 = _result_for(name)
    mismatches = []
    for g in _qpgraph_reference(name):
        # the f3 statistics do not depend on the admixture weights
        r = pm.qpgraph(f2, [tuple(e) for e in g["edges"]], numstart=1)
        if r.pops != g["pops"]:
            mismatches.append(f"{g['name']}: pops {r.pops} != {g['pops']}")
            continue
        for x, row in zip(r.f3, g["f3"]):
            where = f"{g['name']} f3({row['pop1']}; {row['pop2']}, {row['pop3']})"
            if (x.pop1, x.pop2, x.pop3) != (row["pop1"], row["pop2"], row["pop3"]):
                mismatches.append(f"{where}: popmix has f3({x.pop1}; {x.pop2}, {x.pop3})")
            if not np.isclose(x.est, row["est"], rtol=F_RTOL, atol=F_ATOL):
                mismatches.append(f"{where}: est {x.est!r} != {row['est']!r}")
            if not np.isclose(x.se, row["se"], rtol=SE_RTOL, atol=0):
                mismatches.append(f"{where}: se {x.se!r} != {row['se']!r}")
        ppinv, ref = np.array(r.ppinv), np.array(g["ppinv"])
        if not np.allclose(ppinv, ref, rtol=SE_RTOL, atol=SE_RTOL * np.abs(ref).max()):
            mismatches.append(f"{g['name']}: ppinv differs by {np.abs(ppinv - ref).max():.3g}")
    assert not mismatches, "\n".join(mismatches)


def _check_qpgraph_fixed(name: str):
    """At fixed admixture weights a qpgraph fit has no optimizer in it: edge
    lengths, fitted f3 and score are exact functions of the f3 statistics."""
    _, f2 = _result_for(name)
    mismatches = []
    for g in _qpgraph_reference(name):
        edges = [tuple(e) for e in g["edges"]]
        for fx in g["fixed"]:
            w = _as_list(fx["admix_weights"])
            r = pm.qpgraph(f2, edges, admix_weights=w)
            where = f"{g['name']} at weights {w}"
            ours = np.array([e.weight for e in r.edges])
            ref = np.array(fx["weight"], dtype=float)
            # lengths at the zero bound: quadprog leaves rounding noise there
            if not np.allclose(ours, ref, rtol=F_RTOL, atol=F_RTOL * np.abs(ref).max()):
                worst = int(np.argmax(np.abs(ours - ref)))
                mismatches.append(f"{where}: edge {edges[worst]} {ours[worst]!r} != {ref[worst]!r}")
            fit = np.array([f.fit for f in r.f3])
            if not np.allclose(fit, fx["f3_fit"], rtol=F_RTOL, atol=F_ATOL):
                mismatches.append(f"{where}: f3 fit differs by {np.abs(fit - fx['f3_fit']).max():.3g}")
            if not np.isclose(r.score, fx["score"], rtol=F_RTOL, atol=0):
                mismatches.append(f"{where}: score {r.score!r} != {fx['score']!r}")
    assert not mismatches, "\n".join(mismatches)


def _check_qpgraph_fit(name: str):
    """qpgraph() with its admixture weights optimized.

    popmix draws admixtools2's random starts (R's set.seed and runif) and
    runs R's L-BFGS-B from each, so every start must match: the start
    exactly, its final score to the f-statistic tolerance and its weights to
    the weight tolerance (the score differs from admixtools2's by up to 1e-11
    in rounding, which moves the finite-difference path that much). The best
    fit and the low/high range of the weights over the starts follow.

    A graph where either tool's line search failed from some start
    (convergence 52) is compared more loosely. In ancient_5pop admix1 the
    edge s -> Sardinian.DG reaches its zero bound at weight 0.4508, which
    puts a kink in the score with its minimum at 0.45087 (0.677059). The
    finite-difference gradient straddles the kink, so from every start both
    tools stop short of it, at weights 0.4502-0.4503 and scores 0.677125-
    0.677129, at a point set by rounding at the 1e-11 level: the scores then
    agree to STALL_RTOL, the weights to STALL_WEIGHT_TOL and the
    convergence codes not at all.
    """
    _, f2 = _result_for(name)
    mismatches = []
    for g in _qpgraph_reference(name):
        r = pm.qpgraph(f2, [tuple(e) for e in g["edges"]])
        stalled = any(o.convergence == 52 for o in r.opt) or any(row["convergence"] == 52 for row in g["opt"] or [])
        score_rtol, weight_tol = (STALL_RTOL, STALL_WEIGHT_TOL) if stalled else (F_RTOL, WEIGHT_TOL)
        if not np.isclose(r.score, g["score"], rtol=score_rtol, atol=0):
            mismatches.append(f"{g['name']}: score {r.score!r} != admixtools2's {g['score']!r}")
        assert len(r.opt) == len(g["opt"] or [])
        for k, (o, row) in enumerate(zip(r.opt, g["opt"] or [])):  # noqa: B905
            where = f"{g['name']} start {k}"
            if o.start != [row[c] for c in row if c.startswith("i.")]:
                mismatches.append(f"{where}: start {o.start} differs")
            if not np.isclose(o.score, row["value"], rtol=score_rtol, atol=0):
                mismatches.append(f"{where}: score {o.score!r} != {row['value']!r}")
            if not np.allclose(o.weights, [row[c] for c in row if c.startswith("e.")], rtol=0, atol=weight_tol):
                mismatches.append(f"{where}: weights {o.weights} differ")
            if not stalled and o.convergence != row["convergence"]:
                mismatches.append(f"{where}: convergence {o.convergence} != {row['convergence']}")
        ref_e = {(e["from"], e["to"]): e for e in g["fit_edges"] if e["type"] == "admix"}
        for e in r.edges:
            if e.kind != "admix":
                continue
            theirs = ref_e[(e.from_node, e.to_node)]
            for field in ("weight", "low", "high"):
                if abs(getattr(e, field) - theirs[field]) > weight_tol:
                    mismatches.append(f"{g['name']} {e.from_node}->{e.to_node}: {field} "
                                      f"{getattr(e, field)!r} != {theirs[field]!r}")
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
            pseudohaploid="admixtools2",
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
    "qpadm_popdrop_corrected": _check_qpadm_popdrop_corrected,
    "f4ratio": _check_f4ratio,
    "qpgraph_f3": _check_qpgraph_f3,
    "qpgraph_fixed": _check_qpgraph_fixed,
    "qpgraph_fit": _check_qpgraph_fit,
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
