"""Many qpAdm models at once, as admixtools2's ``qpadm_p``, ``test_cladality``,
``qpwave_pairs``, ``qpadm_rotate``, ``rotate_models``, ``qpadm_multi`` and
``qpadm_models``.

Each takes and returns what its admixtools2 namesake does: models are
generated in admixtools2's order, and tables are data frames with one model
per row (populations in list columns).
"""

from __future__ import annotations

import itertools
from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd

from popmix import _core


def _pops(pops: str | Sequence[str] | None) -> list[str]:
    if pops is None:
        return []
    return [pops] if isinstance(pops, str) else list(pops)


def _power_set(pops: Sequence[str]) -> list[list[str]]:
    """admixtools2 ``power_set``: subsets by size, then in ``combn`` order."""
    return [list(c) for k in range(1, len(pops) + 1) for c in itertools.combinations(pops, k)]


def qpadm_p(
    f2_blocks: _core.PyF2Blocks,
    left: str | Sequence[str],
    right: str | Sequence[str],
    target: str | None = None,
    *,
    rnk: int | None = None,
    fudge: float = 1e-4,
    constrained: bool = False,
    weights: bool = False,
) -> dict[str, Any]:
    """Fit one qpAdm model at one rank, without standard errors or popdrop.

    Parameters
    ----------
    f2_blocks
        f2 blocks with every population.
    left, right
        Left and right populations. With a ``target`` it is put first among
        the left populations; otherwise the first left population is the
        base.
    target
        Target population.
    rnk
        Rank of the fit (default: the number of left populations other than
        the target, minus 1: the qpAdm model).
    fudge, constrained
        As for ``qpadm``.
    weights
        Also return the weights.

    Returns
    -------
    dict
        ``f4rank``, ``dof``, ``chisq``, ``p`` (None when ``dof`` is 0),
        ``feasible`` (all weights in [0, 1]) and, with ``weights=True``,
        ``weights``: population to weight, for the left populations after the
        first, or None at ranks where the fit does not determine them.
    """
    left, right = _pops(left), _pops(right)
    r = _core.qpadm_p(f2_blocks, left, right, target, rnk, fudge, constrained)
    w = r.pop("weights")
    if weights:
        rows = ([target] if target is not None else []) + [p for p in left if p != target]
        r["weights"] = None if w is None else dict(zip(rows[1:], w))  # noqa: B905
    return r


def test_cladality(
    f2_blocks: _core.PyF2Blocks, left: Sequence[str], right: Sequence[str], fudge: float = 1e-4
) -> dict[str, Any]:
    """Whether ``left`` forms a clade with respect to ``right``: ``qpadm_p`` at rank 0."""
    return qpadm_p(f2_blocks, left, right, rnk=0, fudge=fudge)


# not a pytest test, whatever pytest makes of the name
test_cladality.__test__ = False  # type: ignore[attr-defined]


def qpwave_pairs(
    f2_blocks: _core.PyF2Blocks, left: Sequence[str], right: Sequence[str], fudge: float = 1e-4
) -> pd.DataFrame:
    """Rank-0 qpWave chi-squared and p-value for every pair of ``left``
    populations: whether each pair forms a clade with respect to ``right``.

    Returns
    -------
    pandas.DataFrame
        ``pop1``, ``pop2``, ``chisq``, ``p``; each pair in both orders,
        sorted by ``pop1`` then ``pop2``.
    """
    rows = []
    for a, b in itertools.combinations(sorted(set(left)), 2):
        r = qpadm_p(f2_blocks, [b], right, a, rnk=0, fudge=fudge)
        rows += [(a, b, r["chisq"], r["p"]), (b, a, r["chisq"], r["p"])]
    return pd.DataFrame(sorted(rows), columns=["pop1", "pop2", "chisq", "p"])


def _all_lr2(pops: Sequence[str], rightfix: int) -> list[tuple[list[str], list[str]]]:
    """admixtools2 ``all_lr2``: every left set of up to (n + rightfix - 1) / 2
    populations, and the rest on the right."""
    maxleft = min(len(pops), (len(pops) + rightfix - 1) // 2)
    lefts = [list(c) for k in range(1, maxleft + 1) for c in itertools.combinations(pops, k)]
    return [(left, [p for p in pops if p not in left]) for left in lefts]


def rotate_models(
    leftright: Sequence[str], target: str, rightfix: str | Sequence[str] | None = None
) -> pd.DataFrame:
    """Models that distribute ``leftright`` between left and right, as
    admixtools2 ``rotate_models``: the input for ``qpadm_multi``.

    Returns
    -------
    pandas.DataFrame
        ``left``, ``right`` (the populations not on the left, then
        ``rightfix``) and ``target``.
    """
    fix = _pops(rightfix)
    return pd.DataFrame(
        [(left, right + fix, target) for left, right in _all_lr2(list(leftright), len(fix))],
        columns=["left", "right", "target"],
    )


def qpadm_rotate(
    f2_blocks: _core.PyF2Blocks,
    leftright: Sequence[str],
    target: str,
    rightfix: str | Sequence[str] | None = None,
    full_results: bool = False,
    fudge: float = 1e-4,
    constrained: bool = False,
) -> pd.DataFrame | list[_core.PyQpAdmResult]:
    """Every qpAdm model that keeps some of ``leftright`` as sources of
    ``target`` and moves the rest to the right, next to ``rightfix``
    (admixtools2 ``qpadm_rotate``).

    Returns
    -------
    pandas.DataFrame or list of QpAdmResult
        One row per model: ``left``, ``right`` (``rightfix`` first), and the
        ``qpadm_p`` columns; with ``full_results``, each model's full
        ``qpadm`` result instead.
    """
    fix = _pops(rightfix)
    models = [(left, fix + right) for left, right in _all_lr2(list(leftright), len(fix))]
    if full_results:
        return [_core.qpadm(f2_blocks, target, left, right, constrained, fudge) for left, right in models]
    rows = [{"left": left, "right": right,
             **qpadm_p(f2_blocks, left, right, target, fudge=fudge, constrained=constrained)}
            for left, right in models]
    return pd.DataFrame(rows, columns=["left", "right", "f4rank", "dof", "chisq", "p", "feasible"])


def _model_rows(models: pd.DataFrame | Mapping[str, Sequence] | Sequence[Mapping[str, Any]]) -> list[dict]:
    if isinstance(models, pd.DataFrame):
        records = models.to_dict("records")
    elif isinstance(models, Mapping):
        n = {len(v) for v in models.values()}
        if len(n) != 1:
            raise ValueError("'left', 'right' (and 'target') should be of equal length")
        records = [dict(zip(models, row)) for row in zip(*models.values())]  # noqa: B905
    else:
        records = [dict(m) for m in models]
    out = []
    for m in records:
        if "left" not in m or "right" not in m:
            raise ValueError("each model needs 'left' and 'right' populations")
        out.append({"left": _pops(m["left"]), "right": _pops(m["right"]), "target": m.get("target")})
    return out


def qpadm_multi(
    f2_blocks: _core.PyF2Blocks,
    models: pd.DataFrame | Mapping[str, Sequence] | Sequence[Mapping[str, Any]],
    full_results: bool = True,
    fudge: float = 1e-4,
    constrained: bool = False,
) -> pd.DataFrame | list[_core.PyQpAdmResult]:
    """Run several qpAdm models (admixtools2 ``qpadm_multi``).

    Parameters
    ----------
    models
        Models with ``left``, ``right`` and optionally ``target``
        populations: a data frame (one model per row, such as
        ``rotate_models`` returns), a dict of equal-length lists, or a list
        of dicts. Without a target, the first left population is the target.

    Returns
    -------
    list of QpAdmResult or pandas.DataFrame
        Each model's ``qpadm`` result; with ``full_results=False``, a data
        frame with the models (``left`` without the target) and the
        ``qpadm_p`` columns.
    """
    rows = []
    for m in _model_rows(models):
        target = m["target"]
        left = [p for p in m["left"] if p != target]
        if target is None:
            target, left = left[0], left[1:]
        if full_results:
            rows.append(_core.qpadm(f2_blocks, target, left, m["right"], constrained, fudge))
        else:
            rows.append({"left": left, "right": m["right"], "target": m["target"],
                         **qpadm_p(f2_blocks, left, m["right"], target, fudge=fudge, constrained=constrained)})
    if full_results:
        return rows
    return pd.DataFrame(rows, columns=["left", "right", "target", "f4rank", "dof", "chisq", "p", "feasible"])


def qpadm_models(pops: Sequence[str], allpops: bool = True, more_right: bool = True) -> pd.DataFrame:
    """Every split of ``pops`` into left and right populations (admixtools2
    ``qpadm_models``).

    Parameters
    ----------
    allpops
        Use every population in each model (the rest of ``pops`` on the
        right); otherwise any non-empty subset of the rest.
    more_right
        Only models with more right than left populations.

    Returns
    -------
    pandas.DataFrame
        ``left`` and ``right``, one model per row.
    """
    pops = list(pops)
    rows = []
    for left in _power_set(pops):
        rest = [p for p in pops if p not in left]
        for right in ([rest] if allpops else _power_set(rest)):
            if right and (not more_right or len(right) > len(left)):
                rows.append((left, right))
    return pd.DataFrame(rows, columns=["left", "right"])
