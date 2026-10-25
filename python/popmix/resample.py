"""Models evaluated on resampled SNP blocks, as admixtools2's
``*_resample_snps``, ``qpgraph_resample_multi``, ``compare_fits``,
``loo_list`` and ``boo_list``.

Resampling leaves out one block at a time (jackknife) or draws blocks with
replacement (bootstrap). Bootstrap draws use R's random number generator:
with ``boot_seed=s`` they are the blocks admixtools2 draws after
``set.seed(s)``. Block indices (``sel``) are 0-based; admixtools2's are
1-based.

As in admixtools2, a model that fails on a resample is recorded in that
row's ``error`` column instead of stopping the run.
"""

from __future__ import annotations

import math
import secrets
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np
import pandas as pd

from popmix import _core
from popmix.dstats import qp3pop, qpdstat


def loo_list(f2_blocks: _core.PyF2Blocks) -> list[_core.PyF2Blocks]:
    """The f2 blocks with each block left out in turn (admixtools2 ``loo_list``)."""
    n = f2_blocks.n_blocks
    return [f2_blocks.select_blocks([b for b in range(n) if b != i]) for i in range(n)]


def _boot_draws(n: int, nboot: int, seed: int | None) -> list[list[int]]:
    if seed is None:
        seed = secrets.randbelow(2**31 - 1)
    return _core.r_sample_replace(seed, n, n, nboot)


def boo_list(f2_blocks: _core.PyF2Blocks, nboot: int | None = None, seed: int | None = None) -> dict[str, list]:
    """Bootstrap resamples of the f2 blocks (admixtools2 ``boo_list``).

    Parameters
    ----------
    nboot
        Number of resamples (default: the number of blocks).
    seed
        R seed: the same draws as admixtools2 after ``set.seed(seed)``
        (default: a random one).

    Returns
    -------
    dict
        ``boo``: each resample (as many blocks as the data, drawn with
        replacement); ``sel``: the blocks drawn (0-based); ``test``: the
        blocks not drawn.
    """
    n = f2_blocks.n_blocks
    sel = _boot_draws(n, n if nboot is None else nboot, seed)
    test = [f2_blocks.select_blocks([b for b in range(n) if b not in set(s)]) for s in sel]
    return {"boo": [f2_blocks.select_blocks(s) for s in sel], "sel": sel, "test": test}


def _resamples(f2_blocks: _core.PyF2Blocks, boot: bool | int, boot_seed: int | None) -> tuple[list, list]:
    if boot:
        nboot = f2_blocks.n_blocks if boot is True or boot == 1 else int(boot)
        b = boo_list(f2_blocks, nboot, boot_seed)
        return b["boo"], b["sel"]
    n = f2_blocks.n_blocks
    return loo_list(f2_blocks), [[b for b in range(n) if b != i] for i in range(n)]


def _run(fun: Callable[[Any], dict], data: Sequence, sel: Sequence) -> pd.DataFrame:
    rows = []
    for i, (d, s) in enumerate(zip(data, sel), start=1):  # noqa: B905
        try:
            out, error = fun(d), None
        except Exception as e:  # recorded per row, as admixtools2's safely()
            out, error = {}, e
        rows.append({"id": i, "sel": list(s), **out, "error": error})
    return pd.DataFrame(rows)


def _qpgraph_row(r: _core.PyQpGraphResult) -> dict:
    return {"score": r.score, "score_test": r.score_test,
            "weights": [e.weight for e in r.edges if e.kind == "admix"], "result": r}


def qpgraph_resample_snps(
    f2_blocks: _core.PyF2Blocks,
    edges: Sequence[tuple[str, str]],
    boot: bool | int = False,
    boot_seed: int | None = None,
    **qpgraph_args: Any,
) -> pd.DataFrame:
    """Fit a graph on each resample of the blocks (admixtools2
    ``qpgraph_resample_snps``).

    Parameters
    ----------
    boot
        False (default): leave out each block in turn. True: as many
        bootstrap resamples as blocks; an integer: that many.
    boot_seed
        R seed for the bootstrap draws.
    **qpgraph_args
        Passed to ``qpgraph`` (``seed`` makes its random starts admixtools2's
        with the same seed).

    Returns
    -------
    pandas.DataFrame
        One row per resample: ``id``, ``sel`` (blocks used), ``score``,
        ``weights`` (of the admixture edges, in edge order), ``result`` (the
        qpgraph result) and ``error``.
    """
    data, sel = _resamples(f2_blocks, boot, boot_seed)
    return _run(lambda d: _qpgraph_row(_core.qpgraph(d, list(edges), **qpgraph_args)), data, sel)


def qpgraph_resample_multi(
    f2_blocks: _core.PyF2Blocks,
    graphs: Sequence[Sequence[tuple[str, str]]],
    nboot: int,
    boot_seed: int | None = None,
    **qpgraph_args: Any,
) -> list[pd.DataFrame]:
    """Fit several graphs on the same bootstrap resamples, each scored on the
    blocks its resample left out (admixtools2 ``qpgraph_resample_multi``);
    compare two with ``compare_fits(fits[0].score_test, fits[1].score_test)``.

    Returns
    -------
    list of pandas.DataFrame
        One per graph, one row per resample: as ``qpgraph_resample_snps``,
        with ``score_test``.
    """
    b = boo_list(f2_blocks, nboot, boot_seed)

    def fit(edges: Sequence[tuple[str, str]]) -> pd.DataFrame:
        it = iter(b["test"])
        return _run(lambda d: _qpgraph_row(_core.qpgraph(d, list(edges), f2_blocks_test=next(it), **qpgraph_args)),
                    b["boo"], b["sel"])

    return [fit(g) for g in graphs]


def compare_fits(scores1: Sequence[float], scores2: Sequence[float]) -> dict[str, float]:
    """Whether one graph fits resamples better than another (admixtools2
    ``compare_fits``), from their scores on the same bootstrap resamples.

    Pairs with a missing or zero difference are left out.

    Returns
    -------
    dict
        ``diff`` (mean score difference), ``se``, ``z``, ``p`` (normal),
        ``p_emp`` (two-sided bootstrap p-value, at least 2 / (n + 1)),
        ``p_emp_nocorr`` (without that floor), ``ci_low`` and ``ci_high``
        (2.5% and 97.5% quantiles of the differences).
    """
    if len(scores1) != len(scores2):
        raise ValueError("Length of scores1 and scores2 differ!")
    d = np.asarray(scores1, dtype=float) - np.asarray(scores2, dtype=float)
    d = d[~np.isnan(d) & (d != 0)]
    n = len(d)
    diff = d.mean()
    # admixtools2 boot_mat_stats with unit weights: the sample variance
    se = math.sqrt(((d - diff) ** 2).sum() / (n - 1))
    z = diff / se
    frac1, frac2 = np.mean(d < 0), np.mean(d > 0)
    return {
        "diff": diff,
        "se": se,
        "z": z,
        "p": math.erfc(abs(z) / math.sqrt(2)),  # 2 pnorm(-|z|)
        "p_emp": min((np.sum(d > 0) + 1) / (n + 1), (np.sum(d < 0) + 1) / (n + 1)) * 2,
        "p_emp_nocorr": min(frac1, frac2) * 2,
        "ci_low": float(np.quantile(d, 0.025)),
        "ci_high": float(np.quantile(d, 0.975)),
    }


def qpadm_resample_snps(
    f2_blocks: _core.PyF2Blocks,
    target: str,
    left: Sequence[str],
    right: Sequence[str],
    boot: bool | int = False,
    boot_seed: int | None = None,
    **qpadm_args: Any,
) -> pd.DataFrame:
    """qpAdm on each resample of the blocks (admixtools2 ``qpadm_resample_snps``).

    Returns
    -------
    pandas.DataFrame
        One row per resample: ``id``, ``sel``, ``weights``, ``chisq`` and
        ``p`` of the full model, ``result`` (the qpadm result) and ``error``.
    """
    def fit(d):
        r = _core.qpadm(d, target, list(left), list(right), **qpadm_args)
        return {"weights": list(r.weights), "chisq": r.chi_sq, "p": r.p_value, "result": r}

    data, sel = _resamples(f2_blocks, boot, boot_seed)
    return _run(fit, data, sel)


def qpwave_resample_snps(
    f2_blocks: _core.PyF2Blocks,
    left: Sequence[str],
    right: Sequence[str],
    boot: bool | int = False,
    boot_seed: int | None = None,
    **qpwave_args: Any,
) -> pd.DataFrame:
    """qpWave on each resample of the blocks (admixtools2 ``qpwave_resample_snps``).

    Returns
    -------
    pandas.DataFrame
        One row per resample: ``id``, ``sel``, ``rankdrop`` (one dict per
        rank, highest first as in admixtools2: ``f4rank``, ``dof``,
        ``chisq``, ``p``), ``result`` (the qpwave result) and ``error``.
    """
    def fit(d):
        r = _core.qpwave(d, list(left), list(right), **qpwave_args)
        rankdrop = [{"f4rank": t.rank, "dof": t.df, "chisq": t.chi_sq, "p": t.p_value}
                    for t in sorted(r.rank_tests, key=lambda t: -t.rank) if t.rank < len(left) - 1]
        return {"rankdrop": rankdrop, "result": r}

    data, sel = _resamples(f2_blocks, boot, boot_seed)
    return _run(fit, data, sel)


def _fstat_resample(fun: Callable, f2_blocks, boot, boot_seed, pops: dict, args: dict) -> pd.DataFrame:
    data, sel = _resamples(f2_blocks, boot, boot_seed)
    frames = []
    for i, (d, s) in enumerate(zip(data, sel), start=1):  # noqa: B905
        t = fun(d, **pops, **args)
        frames.append(t.assign(id=i, sel=[list(s)] * len(t)))
    out = pd.concat(frames, ignore_index=True)
    return out[["id", "sel"] + [c for c in out.columns if c not in ("id", "sel")]]


def qpdstat_resample_snps(
    f2_blocks: _core.PyF2Blocks,
    pop1: Any = None,
    pop2: Any = None,
    pop3: Any = None,
    pop4: Any = None,
    boot: bool | int = False,
    boot_seed: int | None = None,
    **qpdstat_args: Any,
) -> pd.DataFrame:
    """f4-statistics on each resample of the blocks (admixtools2
    ``qpdstat_resample_snps``): ``qpdstat``'s rows for each resample, with its
    ``id`` and ``sel``."""
    return _fstat_resample(qpdstat, f2_blocks, boot, boot_seed,
                           {"pop1": pop1, "pop2": pop2, "pop3": pop3, "pop4": pop4}, qpdstat_args)


def qp3pop_resample_snps(
    f2_blocks: _core.PyF2Blocks,
    pop1: Any = None,
    pop2: Any = None,
    pop3: Any = None,
    boot: bool | int = False,
    boot_seed: int | None = None,
    **qp3pop_args: Any,
) -> pd.DataFrame:
    """f3-statistics on each resample of the blocks (admixtools2
    ``qp3pop_resample_snps``): ``qp3pop``'s rows for each resample, with its
    ``id`` and ``sel``."""
    return _fstat_resample(qp3pop, f2_blocks, boot, boot_seed, {"pop1": pop1, "pop2": pop2, "pop3": pop3},
                           qp3pop_args)
