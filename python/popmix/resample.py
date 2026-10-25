"""Models evaluated on resampled SNP blocks or with individuals left out, as
admixtools2's ``*_resample_snps``, ``*_resample_inds``,
``qpgraph_resample_multi``, ``compare_fits``, ``loo_list`` and ``boo_list``.

Resampling leaves out one block at a time (jackknife) or draws blocks with
replacement (bootstrap). Bootstrap draws use R's random number generator:
with ``boot_seed=s`` they are the blocks admixtools2 draws after
``set.seed(s)``. Block indices (``sel``) are 0-based; admixtools2's are
1-based.

Leaving an individual out recomputes f2 from the genotypes without it, the
same way as ``extract_f2`` (admixtools2's ``*_resample_inds`` read counts
precomputed by ``extract_counts``, a different f2 estimator).

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


def _qpadm_fit(target, left, right, args):
    def fit(d):
        r = _core.qpadm(d, target, list(left), list(right), **args)
        return {"weights": list(r.weights), "chisq": r.chi_sq, "p": r.p_value, "result": r}
    return fit


def _qpwave_fit(left, right, args):
    def fit(d):
        r = _core.qpwave(d, list(left), list(right), **args)
        rankdrop = [{"f4rank": t.rank, "dof": t.df, "chisq": t.chi_sq, "p": t.p_value}
                    for t in sorted(r.rank_tests, key=lambda t: -t.rank) if t.rank < len(left) - 1]
        return {"rankdrop": rankdrop, "result": r}
    return fit


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
    """Compare two graphs' scores on the same bootstrap replicates.

    The bootstrap test of Maier et al. 2023 (eLife, "comparing graphs"),
    with the outputs documented for admixtools2's ``compare_fits``; see also
    ``docs/cleanroom/small-helpers.md`` section 4.

    The differences ``d = scores1 - scores2`` are taken per replicate;
    missing (NaN) differences and exact ties (``d == 0``) are dropped,
    infinite ones kept. With ``n`` differences left, ``n_pos`` of them
    positive and ``n_neg`` negative:

    - ``diff``: the mean of ``d``; ``se``: its standard deviation
      (denominator ``n - 1``); ``z = diff / se``; ``p``: two-sided normal p
      of ``z``.
    - ``p_emp_nocorr``: ``2 * min(n_pos, n_neg) / n``.
    - ``p_emp``: ``2 * min(n_pos + 1, n_neg + 1) / (n + 1)``, so never 0.
    - ``ci_low``, ``ci_high``: the 2.5% and 97.5% quantiles of ``d``
      (linear interpolation between order statistics).

    Parameters
    ----------
    scores1, scores2 : sequence of float
        Scores of two models on the same bootstrap replicates, in the same
        order.

    Returns
    -------
    dict
        ``diff``, ``se``, ``z``, ``p``, ``p_emp``, ``p_emp_nocorr``,
        ``ci_low``, ``ci_high``.

    Raises
    ------
    ValueError
        If the two score vectors differ in length.
    """
    a = np.asarray(scores1, dtype=float)
    b = np.asarray(scores2, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"score vectors differ in length: {a.size} and {b.size}")
    d = a - b
    d = d[~np.isnan(d) & (d != 0)]
    n = d.size
    n_pos = int(np.count_nonzero(d > 0))
    n_neg = int(np.count_nonzero(d < 0))
    nan = float("nan")
    with np.errstate(invalid="ignore", divide="ignore"):
        diff = float(d.mean()) if n else nan
        se = float(np.sqrt(np.sum((d - diff) ** 2) / (n - 1))) if n > 1 else nan
        z = float(np.float64(diff) / se)
    p = nan if math.isnan(z) else math.erfc(abs(z) / math.sqrt(2.0))
    q = np.quantile(d, [0.025, 0.975]) if n else (nan, nan)
    return {
        "diff": diff,
        "se": se,
        "z": z,
        "p": p,
        "p_emp": 2.0 * min(n_pos + 1, n_neg + 1) / (n + 1),
        "p_emp_nocorr": 2.0 * min(n_pos, n_neg) / n if n else nan,
        "ci_low": float(q[0]),
        "ci_high": float(q[1]),
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
    data, sel = _resamples(f2_blocks, boot, boot_seed)
    return _run(_qpadm_fit(target, left, right, qpadm_args), data, sel)


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
    data, sel = _resamples(f2_blocks, boot, boot_seed)
    return _run(_qpwave_fit(left, right, qpwave_args), data, sel)


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


def _left_out(prefix: str, inds: Sequence[str], pops: Sequence[str], f2_args: dict | None):
    """(individual, its population, f2 blocks without it) for every individual
    in a population of more than one."""
    if len(inds) != len(pops):
        raise ValueError("inds and pops must have the same length: one population per individual")
    order = list(dict.fromkeys(pops))
    for ind, pop in zip(inds, pops):  # noqa: B905
        if list(pops).count(pop) < 2:
            continue
        keep = [(i, p) for i, p in zip(inds, pops) if i != ind]  # noqa: B905
        yield ind, pop, lambda keep=keep: _core.extract_f2(
            prefix, order, inds=[i for i, _ in keep], pops=[p for _, p in keep], **(f2_args or {}))


def _run_inds(fun: Callable[[Any], dict], prefix, inds, pops, f2_args) -> pd.DataFrame:
    rows = []
    for ind, pop, f2 in _left_out(prefix, inds, pops, f2_args):
        try:
            out, error = fun(f2()), None
        except Exception as e:  # recorded per row, as admixtools2's safely()
            out, error = {}, e
        rows.append({"ind": ind, "pop": pop, **out, "error": error})
    return pd.DataFrame(rows)


_INDS_DOC = """

    Parameters
    ----------
    prefix
        Genotype prefix.
    inds, pops
        The individuals to use and the population of each (as
        ``extract_f2``'s ``inds`` and ``pops``); each one in a population
        of more than one is left out in turn.
    f2_args
        Passed to ``extract_f2`` (``auto_only``, ``block_size``, ``maxmiss``,
        ``poly_only``).

    Returns
    -------
    pandas.DataFrame
        One row per individual left out: ``ind``, ``pop``, the model's
        columns as in the ``_resample_snps`` function, and ``error``.
    """


def qpadm_resample_inds(prefix: str, inds: Sequence[str], pops: Sequence[str], target: str,
                        left: Sequence[str], right: Sequence[str], f2_args: dict | None = None,
                        **qpadm_args: Any) -> pd.DataFrame:
    return _run_inds(_qpadm_fit(target, left, right, qpadm_args), prefix, inds, pops, f2_args)


def qpwave_resample_inds(prefix: str, inds: Sequence[str], pops: Sequence[str], left: Sequence[str],
                         right: Sequence[str], f2_args: dict | None = None, **qpwave_args: Any) -> pd.DataFrame:
    return _run_inds(_qpwave_fit(left, right, qpwave_args), prefix, inds, pops, f2_args)


def qpgraph_resample_inds(prefix: str, inds: Sequence[str], pops: Sequence[str],
                          edges: Sequence[tuple[str, str]], f2_args: dict | None = None,
                          **qpgraph_args: Any) -> pd.DataFrame:
    return _run_inds(lambda d: _qpgraph_row(_core.qpgraph(d, list(edges), **qpgraph_args)), prefix, inds, pops,
                     f2_args)


def _fstat_inds(fun, prefix, inds, pops, f2_args, pop_args, args) -> pd.DataFrame:
    frames = []
    for ind, pop, f2 in _left_out(prefix, inds, pops, f2_args):
        t = fun(f2(), **pop_args, **args)
        frames.append(t.assign(ind=ind, pop=pop))
    out = pd.concat(frames, ignore_index=True)
    return out[["ind", "pop"] + [c for c in out.columns if c not in ("ind", "pop")]]


def qpdstat_resample_inds(prefix: str, inds: Sequence[str], pops: Sequence[str], pop1: Any = None,
                          pop2: Any = None, pop3: Any = None, pop4: Any = None, f2_args: dict | None = None,
                          **qpdstat_args: Any) -> pd.DataFrame:
    return _fstat_inds(qpdstat, prefix, inds, pops, f2_args,
                       {"pop1": pop1, "pop2": pop2, "pop3": pop3, "pop4": pop4}, qpdstat_args)


def qp3pop_resample_inds(prefix: str, inds: Sequence[str], pops: Sequence[str], pop1: Any = None,
                         pop2: Any = None, pop3: Any = None, f2_args: dict | None = None,
                         **qp3pop_args: Any) -> pd.DataFrame:
    return _fstat_inds(qp3pop, prefix, inds, pops, f2_args, {"pop1": pop1, "pop2": pop2, "pop3": pop3},
                       qp3pop_args)


for _f, _what in [(qpadm_resample_inds, "qpAdm"), (qpwave_resample_inds, "qpWave"),
                  (qpgraph_resample_inds, "qpGraph"), (qpdstat_resample_inds, "f4-statistics"),
                  (qp3pop_resample_inds, "f3-statistics")]:
    _f.__doc__ = f"{_what} with each individual left out in turn (admixtools2 ``{_f.__name__}``)." + _INDS_DOC
