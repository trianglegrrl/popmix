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
from popmix.fits import qpadm
from popmix.qpadm_assess import warn_multiple_testing
from popmix.results import QpAdmResult


def _pops(pops: str | Sequence[str] | None) -> list[str]:
    if pops is None:
        return []
    return [pops] if isinstance(pops, str) else list(pops)


def _subsets(pops: Sequence[str], max_size: int | None = None) -> list[list[str]]:
    """The non-empty subsets of ``pops`` (up to ``max_size`` members), in the
    order admixtools2's model sets use: smaller subsets first, and subsets of
    one size in lexicographic order of their positions in ``pops``."""
    pops = list(pops)
    top = len(pops) if max_size is None else min(max_size, len(pops))
    out: list[list[str]] = []
    for size in range(1, top + 1):
        out.extend([pops[i] for i in idx] for idx in itertools.combinations(range(len(pops)), size))
    return out


def _rotations(leftright: Sequence[str], rightfix: Sequence[str], fixed_first: bool) -> list[tuple[list, list]]:
    """Left/right splits of ``leftright``, with ``rightfix`` always right.

    Every non-empty left set of at most ``(len(leftright) + len(rightfix) - 1)
    // 2`` populations (so the right set is always larger than the left), in
    ``_subsets`` order. The right set is the rest of ``leftright`` (in its
    order) followed by ``rightfix``, or ``rightfix`` first with
    ``fixed_first``.
    """
    leftright, rightfix = list(leftright), list(rightfix)
    splits = []
    for left in _subsets(leftright, (len(leftright) + len(rightfix) - 1) // 2):
        rest = [p for p in leftright if p not in left]
        splits.append((left, rightfix + rest if fixed_first else rest + rightfix))
    return splits


def qpadm_p(
    f2_blocks: _core.F2Blocks,
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
    f2_blocks: _core.F2Blocks, left: Sequence[str], right: Sequence[str], *, fudge: float = 1e-4
) -> dict[str, Any]:
    """Whether ``left`` forms a clade with respect to ``right``: ``qpadm_p`` at rank 0.

    As admixtools2 ``test_cladality``.

    Parameters
    ----------
    f2_blocks
        f2 blocks with every population.
    left, right
        Left and right populations.
    fudge
        As for ``qpadm``.

    Returns
    -------
    dict
        As ``qpadm_p``: ``f4rank``, ``dof``, ``chisq``, ``p`` and
        ``feasible``; a small ``p`` rejects cladality, and a large one does
        not show it (other right populations can reject it).
    """
    return qpadm_p(f2_blocks, left, right, rnk=0, fudge=fudge)


# not a pytest test, whatever pytest makes of the name
test_cladality.__test__ = False  # type: ignore[attr-defined]


def qpwave_pairs(
    f2_blocks: _core.F2Blocks, left: Sequence[str], right: Sequence[str], *, fudge: float = 1e-4
) -> pd.DataFrame:
    """Test every pair of ``left`` populations for cladality with respect to
    ``right`` (admixtools2 ``qpwave_pairs``): qpWave at rank 0.

    Parameters
    ----------
    f2_blocks
        f2 blocks with every population.
    left
        Populations to test in pairs.
    right
        Right populations.
    fudge
        As for ``qpadm``.

    Returns
    -------
    pandas.DataFrame
        ``pop1``, ``pop2``, ``chisq``, ``p``: each pair in both orders, with
        the same statistics, sorted by ``pop1`` and then ``pop2``.
    """
    left, right = _pops(left), _pops(right)
    rows = []
    for a, b in itertools.combinations(left, 2):
        t = test_cladality(f2_blocks, [a, b], right, fudge=fudge)
        rows.append((a, b, t["chisq"], t["p"]))
        rows.append((b, a, t["chisq"], t["p"]))
    rows.sort(key=lambda r: (r[0], r[1]))
    return pd.DataFrame(rows, columns=["pop1", "pop2", "chisq", "p"])


def rotate_models(
    leftright: Sequence[str], target: str, rightfix: str | Sequence[str] | None = None
) -> pd.DataFrame:
    """The qpAdm models of a rotation (admixtools2 ``rotate_models``), for
    ``qpadm_multi``.

    Each non-empty subset of ``leftright`` with fewer members than the
    resulting right set is a left set; the right set is the rest of
    ``leftright`` followed by ``rightfix``. Smaller left sets come first,
    and left sets of one size in the order of their positions in
    ``leftright`` (the order of admixtools2's ``qpadm_rotate`` details).

    Parameters
    ----------
    leftright
        Populations distributed between left and right.
    target
        Target population, the same in every model.
    rightfix
        Populations on the right in every model.

    Returns
    -------
    pandas.DataFrame
        ``left``, ``right``, ``target``, one model per row.
    """
    rows = [(left, right, target) for left, right in _rotations(_pops(leftright), _pops(rightfix), False)]
    return pd.DataFrame(rows, columns=["left", "right", "target"])


def qpadm_rotate(
    f2_blocks: _core.F2Blocks,
    leftright: Sequence[str],
    target: str,
    rightfix: str | Sequence[str] | None = None,
    *, full_results: bool = False,
    fudge: float = 1e-4,
    constrained: bool = False,
) -> pd.DataFrame | list[QpAdmResult]:
    """Fit every qpAdm model of a rotation (admixtools2 ``qpadm_rotate``).

    The models are those of ``rotate_models``, except that the right set is
    ``rightfix`` followed by the rest of ``leftright``. The p-values are not
    corrected for testing them all: a ``PopmixWarning`` says how many models
    were tested, which a report should give with any model chosen (Harney
    et al. 2021); the table's ``attrs["n_models"]`` holds it too.

    Parameters
    ----------
    leftright
        Populations distributed between left and right.
    target
        Target population.
    rightfix
        Populations on the right in every model.
    full_results
        Return each model's full ``qpadm`` result instead of a table.
    fudge, constrained
        As for ``qpadm``.

    Returns
    -------
    pandas.DataFrame or list of QpAdmResult
        A data frame with ``left``, ``right`` and the ``qpadm_p`` columns
        (``f4rank``, ``dof``, ``chisq``, ``p``, ``feasible``), one model
        per row; with ``full_results=True``, the ``qpadm`` results.
    """
    models = [{"left": left, "right": right, "target": target}
              for left, right in _rotations(_pops(leftright), _pops(rightfix), True)]
    out = _multi(f2_blocks, models, full_results, fudge, constrained)
    warn_multiple_testing(len(out))
    if full_results:
        return out
    table = out.drop(columns="target")
    table.attrs["n_models"] = len(table)
    return table


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
        target = m.get("target")
        # a DataFrame holds a missing target as NaN (pandas 3) or None
        if target is not None and not isinstance(target, str) and pd.isna(target):
            target = None
        out.append({"left": _pops(m["left"]), "right": _pops(m["right"]), "target": target})
    return out


def qpadm_multi(
    f2_blocks: _core.F2Blocks,
    models: pd.DataFrame | Mapping[str, Sequence] | Sequence[Mapping[str, Any]],
    *, full_results: bool = True,
    fudge: float = 1e-4,
    constrained: bool = False,
) -> pd.DataFrame | list[QpAdmResult]:
    """Run several qpAdm models (admixtools2 ``qpadm_multi``).

    As for ``qpadm_rotate``, a ``PopmixWarning`` gives the number of models
    tested (also the table's ``attrs["n_models"]``): their p-values are not
    corrected for multiple testing.

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
    out = _multi(f2_blocks, models, full_results, fudge, constrained)
    warn_multiple_testing(len(out))
    if not full_results:
        out.attrs["n_models"] = len(out)
    return out


def _multi(f2_blocks: _core.F2Blocks, models: Any, full_results: bool, fudge: float,
           constrained: bool) -> pd.DataFrame | list[QpAdmResult]:
    rows = []
    for m in _model_rows(models):
        target = m["target"]
        left = [p for p in m["left"] if p != target]
        if target is None:
            target, left = left[0], left[1:]
        if full_results:
            rows.append(qpadm(f2_blocks, left, m["right"], target, constrained=constrained, fudge=fudge))
        else:
            rows.append({"left": left, "right": m["right"], "target": m["target"],
                         **qpadm_p(f2_blocks, left, m["right"], target, fudge=fudge, constrained=constrained)})
    if full_results:
        return rows
    return pd.DataFrame(rows, columns=["left", "right", "target", "f4rank", "dof", "chisq", "p", "feasible"])


def qpadm_models(pops: Sequence[str], *, allpops: bool = True, more_right: bool = True) -> pd.DataFrame:
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
    for left in _subsets(pops):
        rest = [p for p in pops if p not in left]
        for right in ([rest] if allpops else _subsets(rest)):
            if right and (not more_right or len(right) > len(left)):
                rows.append((left, right))
    return pd.DataFrame(rows, columns=["left", "right"])
