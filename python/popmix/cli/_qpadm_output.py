"""qpAdm's output tables (TASKS.md H4): weights with SEs and z, the rank
tests (rankdrop), popdrop, and the whole result as JSON. Tables are columns
for ``write_columns``, which writes floats at 12 significant digits."""

from __future__ import annotations

import dataclasses
import math

import popmix as pm

from ._common import _json_number


def p_text(p: float) -> str:
    """A p-value for text: 3 significant digits, so a small one shows as
    2.25e-13 rather than 0.0000; one below double precision as <1e-300; NA
    for a saturated model."""
    if math.isnan(p):
        return "NA"
    return "<1e-300" if p == 0 else f"{p:.3g}"


def weight_columns(r: pm.QpAdmResult) -> dict:
    n = len(r.sources)
    return {"target": [r.target] * n, "source": list(r.sources), "weight": list(r.weights),
            "se": list(r.weight_se), "z": list(r.weight_z), "chisq": [r.chisq] * n, "dof": [r.dof] * n,
            "p": [r.p] * n, "feasible": [r.feasible] * n, "n_blocks": [r.n_blocks] * n,
            "n_snps": [r.n_snps] * n}


def rankdrop_columns(r: pm.QpAdmResult) -> dict:
    return {"f4rank": [t.rank for t in r.rankdrop], "dof": [t.dof for t in r.rankdrop],
            "chisq": [t.chisq for t in r.rankdrop], "p": [t.p for t in r.rankdrop]}


def popdrop_columns(r: pm.QpAdmResult) -> dict:
    rows = r.popdrop
    nan = float("nan")
    cols: dict = {"pattern": [x.pattern for x in rows], "wt": [x.wt for x in rows], "dof": [x.dof for x in rows],
                  "chisq": [x.chisq for x in rows], "p": [x.p for x in rows], "f4rank": [x.f4rank for x in rows],
                  "feasible": [x.feasible for x in rows], "best": ["" if x.best is None else x.best for x in rows],
                  "dofdiff": ["" if x.dofdiff is None else x.dofdiff for x in rows],
                  "chisqdiff": [nan if x.chisqdiff is None else x.chisqdiff for x in rows],
                  "p_nested": [nan if x.p_nested is None else x.p_nested for x in rows]}
    for k, source in enumerate(r.sources):
        cols[source] = [nan if x.weights[k] is None else x.weights[k] for x in rows]
    return cols


def _popdrop_json(row: pm.PopDropRow) -> dict:
    out = dataclasses.asdict(row)
    for key in ("chisq", "p", "chisqdiff", "p_nested"):
        out[key] = _json_number(out[key]) if out[key] is not None else None
    out["weights"] = [None if w is None else _json_number(w) for w in row.weights]
    return out


def as_json(r: pm.QpAdmResult, assessment: pm.QpAdmAssessment) -> dict:
    return {
        "target": r.target,
        "sources": list(r.sources),
        "weights": dict(zip(r.sources, r.weights)),
        "se": dict(zip(r.sources, r.weight_se)),
        "z": {s: _json_number(z) for s, z in zip(r.sources, r.weight_z)},
        "chisq": r.chisq,
        "dof": r.dof,
        "p": _json_number(r.p),
        "feasible": r.feasible,
        "n_blocks": r.n_blocks,
        "n_snps": r.n_snps,
        "assessment": dataclasses.asdict(assessment) | {"notes": list(assessment.notes)},
        "rankdrop": [{"f4rank": t.rank, "dof": t.dof, "chisq": t.chisq, "p": _json_number(t.p)} for t in r.rankdrop],
        "popdrop": [_popdrop_json(row) for row in r.popdrop],
    }
