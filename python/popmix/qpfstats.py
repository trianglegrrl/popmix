"""admixtools2's ``qpfstats`` and ``qpf4diff``.

``qpfstats`` estimates f2 blocks that agree with every f2-, f3- and
f4-statistic of a set of populations, each computed from the SNPs with data
in its own populations (admixtools2 ``allsnps = TRUE``): it regresses the
statistics, block by block, on the pairs' f2 (a ridge of 1e-5), and shifts
each pair's blocks so their mean is the regression of the statistics over
all blocks. ``qpf4diff`` is the jackknife of the difference
f4(pop1, pop2; pop3, pop4) - f4(pop1, pop2; pop5, pop4).
"""

from __future__ import annotations

import os
from collections.abc import Sequence

import numpy as np
import pandas as pd

from popmix import _core
from popmix.dstats import _index, fstat_combinations


def _pair_index(pops: Sequence[str]) -> dict[tuple[str, str], int]:
    out, k = {}, 0
    for i in range(len(pops)):
        for j in range(i + 1, len(pops)):
            out[pops[i], pops[j]] = out[pops[j], pops[i]] = k
            k += 1
    return out


def _fstat_matrix(combs: list[tuple], pops: list[str]) -> np.ndarray:
    """admixtools2 ``construct_fstat_matrix``: each statistic as a
    combination of the pairs' f2."""
    idx = _pair_index(pops)
    x = np.zeros((len(combs), len(idx) // 2))
    for r, (a, b, c, d) in enumerate(combs):
        # assignments, not sums, and a population with itself is skipped
        for (p, q), v in (((a, d), 1), ((b, c), 1), ((a, c), -1), ((b, d), -1)):
            if p != q:
                x[r, idx[p, q]] = v
        if a == c and b == d:
            x[r] *= 2
    return x / 2


def qpfstats(prefix: str | os.PathLike, pops: Sequence[str], include_f2: bool = True, include_f3: bool = True,
             include_f4: bool = True, block_size: float = 0.05, auto_only: bool = True) -> _core.PyF2Blocks:
    """f2 blocks fitted to the f2-, f3- and f4-statistics of ``pops``
    (admixtools2 ``qpfstats``).

    Returns
    -------
    F2Blocks
        For the populations in sorted order, on the blocks of all SNPs.
    """
    sp = sorted(pops)
    combs: list[tuple] = []
    if include_f2:
        combs += [(a, b, a, b) for a in sp for b in sp if a < b]
    if include_f3:
        combs += [(a, b, a, c) for a, b, c in fstat_combinations(list(pops), fnum=3)]
    if include_f4:
        combs += list(fstat_combinations(list(pops), fnum=4))
    upops, quads = _index(combs)
    out = _core.f4_from_geno(os.fspath(prefix), upops, quads, True, block_size, auto_only)
    nb = len(out["block_lengths"])
    ymat = np.asarray(out["block_num"], dtype=float).reshape(len(combs), nb)
    ymat[np.asarray(out["block_n"]).reshape(len(combs), nb) == 0] = np.nan
    y = np.asarray(out["est"], dtype=float)
    plist = list(dict.fromkeys(p for col in zip(*combs) for p in col))  # noqa: B905
    x = _fstat_matrix(combs, plist)
    lh = np.linalg.solve(x.T @ x + np.eye(x.shape[1]) * 1e-5, x.T)
    b, bglob = lh @ ymat, lh @ y
    lengths = np.asarray(out["block_lengths"], dtype=float)
    idx = _pair_index(plist)
    n = len(sp)
    means = [[[0.0] * nb for _ in sp] for _ in sp]
    for i in range(n):
        for j in range(n):
            if i != j:
                k = idx[sp[i], sp[j]]
                blocks = b[k]
                # shift the blocks so their mean (over the blocks with data:
                # a block where some statistic has none is NaN) is the
                # all-block fit
                ok = ~np.isnan(blocks)
                shift = bglob[k] - np.sum(blocks[ok] * lengths[ok]) / lengths[ok].sum()
                means[i][j] = list(blocks + shift)
    # blocks without data carry no weight
    counts = [[[0.0 if np.isnan(v) else n for v, n in zip(means[i][j], lengths)] if i != j else list(lengths)  # noqa: B905
               for j in range(n)] for i in range(n)]
    return _core.PyF2Blocks.from_block_means(sp, [int(v) for v in lengths], means, counts)


def _loo(arr: np.ndarray, bl: np.ndarray) -> np.ndarray:
    """admixtools2 ``est_to_loo`` for one statistic's blocks."""
    ok = ~np.isnan(arr)
    tot = np.sum(arr[ok] * bl[ok]) / bl[ok].sum()
    rel = bl / bl[ok].sum()
    return (tot - arr * rel) / (1 - rel)


def _jack(loo: np.ndarray, bl: np.ndarray) -> tuple[float, float]:
    """admixtools2 ``jack_mat_stats`` without a total (``jack_mat_stats2``)."""
    ok = ~np.isnan(loo)
    h = bl[ok].sum() / bl[ok]
    w = 1 - 1 / h
    est = np.sum(loo[ok] * w) / w.sum()
    xtau = (est - loo[ok]) * np.sqrt(h - 1)
    return float(est), float(np.sqrt(np.sum(xtau**2) / ok.sum()))


def qpf4diff(f2_blocks: _core.PyF2Blocks, pops: Sequence[str] | Sequence[Sequence[str]]) -> pd.DataFrame:
    """f4(pop1, pop2; pop3, pop4) - f4(pop1, pop2; pop5, pop4), with its block
    jackknife standard error (admixtools2 ``qpf4diff``), for one set of five
    populations or several.

    Returns
    -------
    pandas.DataFrame
        pop1-pop5, f4diff, se, z and n (SNPs).
    """
    rows = [tuple(pops)] if isinstance(pops[0], str) else [tuple(p) for p in pops]
    if any(len(r) != 5 for r in rows):
        raise ValueError("'pops' should be 5 populations, or a list of 5-population sets")
    bl = np.asarray(f2_blocks.block_lengths, dtype=float)

    def f2(a: str, b: str) -> np.ndarray:
        return np.asarray(f2_blocks.f2_per_block(a, b), dtype=float)

    def f4(a: str, b: str, c: str, d: str) -> np.ndarray:
        return (f2(a, d) + f2(b, c) - f2(a, c) - f2(b, d)) / 2

    out = []
    for p1, p2, p3, p4, p5 in rows:
        est, se = _jack(_loo(f4(p1, p2, p3, p4) - f4(p1, p2, p5, p4), bl), bl)
        out.append((p1, p2, p3, p4, p5, est, se, est / se, int(bl.sum())))
    return pd.DataFrame(out, columns=["pop1", "pop2", "pop3", "pop4", "pop5", "f4diff", "se", "z", "n"])
