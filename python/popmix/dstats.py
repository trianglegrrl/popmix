"""f4-, D- and f3-statistic scans, as admixtools2's ``qpdstat`` and ``qp3pop``.

Given a genotype prefix, statistics are computed from the genotypes, with every
autosomal SNP and each statistic's own SNPs (admixtools2 ``allsnps = TRUE``);
D-statistics (``f4mode=False``) and f3 normalized by the target's
heterozygosity are available only this way. Given f2 blocks, f4 and f3 come
from the blocks (one SNP set for all statistics).
"""

from __future__ import annotations

import itertools
import os
from collections.abc import Sequence
from typing import Union

import numpy as np
import pandas as pd

from popmix import _core

# Union, not X | Y: aliases are evaluated at runtime, and popmix supports Python 3.9
Data = Union[str, "os.PathLike[str]", "_core.PyF2Blocks"]  # noqa: UP007
Pops = Union[str, Sequence[str], None]  # noqa: UP007

def _as_list(pops: Pops) -> list[str] | None:
    if pops is None:
        return None
    return [pops] if isinstance(pops, str) else list(pops)


def _distinct(rows: list[tuple]) -> list[tuple]:
    seen: set = set()
    return [r for r in rows if not (r in seen or seen.add(r))]


def fstat_combinations(
    pop1: Pops,
    pop2: Pops = None,
    pop3: Pops = None,
    pop4: Pops = None,
    *,
    fnum: int,
    unique_only: bool = True,
) -> list[tuple[str, ...]]:
    """The population combinations of an f3 (``fnum=3``) or f4 (``fnum=4``)
    scan, in the order of admixtools2's ``qp3pop`` and ``qpdstat``.

    Established from those functions' documentation (``pop1``: a vector of
    populations, combined with the other ``pop`` arguments; ``unique_only``:
    redundant combinations excluded) and their output; see also
    ``docs/cleanroom/small-helpers.md`` section 4.

    With ``pop1`` alone (repeats dropped), and ``unique_only``:

    - f3: each 3-subset (a, b, c) of ``pop1``, subsets in lexicographic
      order of positions, taken in its three rotations f3(a; b, c),
      f3(b; c, a), f3(c; a, b): all first rotations, then all second, then
      all third;
    - f4: each 4-subset (a, b, c, d) in its three splits (a, b; c, d),
      (a, c; b, d), (a, d; b, c), one subset's splits together.

    Without ``unique_only``, f3 is every ordered triple, repeats included,
    with the first position varying slowest; f4 pairs every 2-subset with
    every 2-subset (both in lexicographic order).

    With every ``pop`` argument given, the combinations are their grid
    (``pop1`` varying slowest), repeated combinations dropped;
    ``unique_only`` does not apply.

    Parameters
    ----------
    pop1, pop2, pop3, pop4
        Population names (a string is one population). ``pop2`` to
        ``pop{fnum}`` are given together or not at all.
    fnum
        3 or 4.
    unique_only
        Leave out redundant combinations (``pop1`` alone only).

    Returns
    -------
    list of tuple of str
        One tuple of ``fnum`` populations per statistic.

    Raises
    ------
    ValueError
        For another ``fnum``, ``pop2``... given only in part, or fewer than
        ``fnum`` distinct populations with ``unique_only``.
    """
    if fnum not in (3, 4):
        raise ValueError(f"fnum must be 3 or 4, not {fnum}")
    others = [pop2, pop3, pop4][: fnum - 1]
    if [p for p in [pop2, pop3, pop4][fnum - 1:] if p is not None]:
        raise ValueError(f"an f{fnum}-statistic takes populations pop1 to pop{fnum} only")
    lists = [_as_list(p) for p in [pop1, *others]]
    if all(p is not None for p in lists[1:]):
        if lists[0] is None:
            raise ValueError("pop1 is needed with the other populations")
        return _distinct(list(itertools.product(*[list(dict.fromkeys(p)) for p in lists])))
    if any(p is not None for p in lists[1:]):
        raise ValueError(f"give pop2 to pop{fnum} together, or pop1 alone")
    pops = list(dict.fromkeys(lists[0] or []))
    if not unique_only:
        if fnum == 3:
            return list(itertools.product(pops, repeat=3))
        pairs = list(itertools.combinations(pops, 2))
        return [(*x, *y) for x in pairs for y in pairs]
    if len(pops) < fnum:
        raise ValueError(f"an f{fnum} scan needs at least {fnum} distinct populations, got {len(pops)}")
    subsets = list(itertools.combinations(pops, fnum))
    if fnum == 3:
        return [(s[r], s[(r + 1) % 3], s[(r + 2) % 3]) for r in range(3) for s in subsets]
    return [split for a, b, c, d in subsets for split in ((a, b, c, d), (a, c, b, d), (a, d, b, c))]


def _is_prefix(data: Data) -> bool:
    return isinstance(data, (str, os.PathLike))


def _index(combs: Sequence[Sequence[str]]) -> tuple[list[str], list[list[int]]]:
    pops = list(dict.fromkeys(p for c in combs for p in c))
    where = {p: i for i, p in enumerate(pops)}
    return pops, [[where[p] for p in c] for c in combs]


def _stats_frame(combs: list[tuple], cols: list[str], out: dict) -> pd.DataFrame:
    df = pd.DataFrame(combs, columns=cols)
    for k in ("est", "se", "z", "p"):
        df[k] = np.asarray(out[k])
    df["n"] = np.asarray(out["n"]).astype(np.int64)
    return df


def _blocks_frame(combs: list[tuple], cols: list[str], out: dict, num: str, den: str) -> pd.DataFrame:
    lengths = np.asarray(out["block_lengths"], dtype=np.int64)
    nb = len(lengths)
    rows = [c for c in combs for _ in range(nb)]
    df = pd.DataFrame(rows, columns=cols)
    df["block"] = np.tile(np.arange(1, nb + 1), len(combs))
    n = np.asarray(out["block_n"]).ravel().astype(np.int64)
    numer = np.asarray(out["block_num"]).ravel()
    denom = np.asarray(out["block_den"]).ravel()
    df[num] = numer
    df[den] = denom
    df["n"] = n
    df["length"] = np.tile(lengths, len(combs))
    return df


def f4blockdat_from_geno(
    prefix: str | os.PathLike[str],
    popcombs: Sequence[Sequence[str]],
    *,
    f4mode: bool = True,
    block_size: float = 0.05,
    auto_only: bool = True,
    pseudohaploid: str = "detect",
) -> pd.DataFrame:
    """Per-block f4 (or D numerator and denominator) from genotypes.

    As admixtools2 ``f4blockdat_from_geno(allsnps = TRUE)``: one row per
    combination and block, with ``est`` (mean numerator over the block's SNPs
    with data in all four populations; NaN if none), ``den`` (mean D
    denominator; NaN in f4 mode), ``n`` (SNPs used) and ``length`` (SNPs in
    the block).
    """
    combs = [tuple(c) for c in popcombs]
    pops, quads = _index(combs)
    out = _core.f4_from_geno(os.fspath(prefix), pops, quads, f4mode, block_size, auto_only, pseudohaploid)
    df = _blocks_frame(combs, ["pop1", "pop2", "pop3", "pop4"], out, "est", "den")
    df.loc[df["n"] == 0, "est"] = np.nan
    return df


def f3blockdat_from_geno(
    prefix: str | os.PathLike[str],
    popcombs: Sequence[Sequence[str]],
    *,
    outgroupmode: bool = False,
    apply_corr: bool = True,
    block_size: float = 0.05,
    auto_only: bool = True,
    pseudohaploid: str = "detect",
) -> pd.DataFrame:
    """Per-block f3 from genotypes, as admixtools2 ``f3blockdat_from_geno(allsnps = TRUE)``.

    One row per combination and block: ``numer`` (mean bias-corrected
    ``(a-b)(a-c)``), ``denom`` (mean twice the target's heterozygosity, or 1
    in outgroup mode), ``est`` (their ratio), ``n`` and ``length``.
    """
    combs = [tuple(c) for c in popcombs]
    pops, triples = _index(combs)
    out = _core.f3_from_geno(os.fspath(prefix), pops, triples, outgroupmode, apply_corr, block_size,
                             auto_only, pseudohaploid)
    df = _blocks_frame(combs, ["pop1", "pop2", "pop3"], out, "numer", "denom")
    df["est"] = df["numer"] / df["denom"]
    df.loc[(df["n"] == 0) | ~np.isfinite(df["est"]), "est"] = np.nan
    return df


def qpdstat(
    data: Data,
    pop1: Pops = None,
    pop2: Pops = None,
    pop3: Pops = None,
    pop4: Pops = None,
    *,
    f4mode: bool = True,
    unique_only: bool = True,
    block_size: float = 0.05,
    auto_only: bool = True,
    pseudohaploid: str = "detect",
) -> pd.DataFrame:
    """f4- or D-statistics f(pop1, pop2; pop3, pop4), as admixtools2 ``qpdstat``.

    Parameters
    ----------
    data
        A genotype prefix, or f2 blocks (``extract_f2``).
    pop1, pop2, pop3, pop4
        Populations (see ``fstat_combinations``). With ``pop1`` alone, every
        distinct statistic on those populations; with f2 blocks and no
        populations, all of the blocks' populations.
    f4mode
        False for D-statistics, ``(a-b)(c-d) / ((a+b-2ab)(c+d-2cd))``; needs a
        genotype prefix (admixtools2 silently returns f4 for f2 blocks; here
        it is an error).
    unique_only
        As ``fstat_combinations``.
    block_size, auto_only
        Jackknife block size (Morgans) and autosomes only; genotypes only.
    pseudohaploid
        From genotypes: which samples count one allele per SNP, as
        ``sample_ploidy``.

    Returns
    -------
    pandas.DataFrame
        pop1..pop4, est, se, z, p and, from genotypes, n (SNPs used).
    """
    if pop1 is None and not _is_prefix(data):
        pop1 = list(data.populations)
    combs = fstat_combinations(pop1, pop2, pop3, pop4, fnum=4, unique_only=unique_only)
    cols = ["pop1", "pop2", "pop3", "pop4"]
    if _is_prefix(data):
        pops, quads = _index(combs)
        out = _core.f4_from_geno(os.fspath(data), pops, quads, f4mode, block_size, auto_only, pseudohaploid)
        return _stats_frame(combs, cols, out)
    if not f4mode:
        raise ValueError("D-statistics (f4mode=False) need genotypes, not f2 blocks")
    return _from_f2(combs, cols, lambda c: data.f4(*c))


def qp3pop(
    data: Data,
    pop1: Pops = None,
    pop2: Pops = None,
    pop3: Pops = None,
    *,
    outgroupmode: bool = False,
    apply_corr: bool = True,
    unique_only: bool = True,
    block_size: float = 0.05,
    auto_only: bool = True,
    pseudohaploid: str = "detect",
) -> pd.DataFrame:
    """f3-statistics f3(pop1; pop2, pop3), as admixtools2 ``qp3pop``.

    Parameters
    ----------
    data
        A genotype prefix, or f2 blocks.
    pop1, pop2, pop3
        Populations (see ``fstat_combinations``); pop1 is the target.
    outgroupmode
        From genotypes: leave f3 unnormalized (outgroup f3). Otherwise, as in
        admixtools2, f3 from genotypes is divided by the target's
        heterozygosity. f3 from f2 blocks is never normalized.
    apply_corr
        From genotypes: subtract the target's sampling bias.
    unique_only, block_size, auto_only, pseudohaploid
        As ``qpdstat``.

    Returns
    -------
    pandas.DataFrame
        pop1..pop3, est, se, z, p and, from genotypes, n.
    """
    if pop1 is None and not _is_prefix(data):
        pop1 = list(data.populations)
    combs = fstat_combinations(pop1, pop2, pop3, fnum=3, unique_only=unique_only)
    return _f3_rows(data, combs, outgroupmode=outgroupmode, apply_corr=apply_corr, block_size=block_size,
                    auto_only=auto_only, pseudohaploid=pseudohaploid)


def _from_f2(combs: list[tuple], cols: list[str], stat) -> pd.DataFrame:
    from statistics import NormalDist

    rows = []
    for c in combs:
        r = stat(c)
        est, se = r["estimate"], r["std_err"]
        z = est / se if se > 0 else float("nan")
        p = 2 * NormalDist().cdf(-abs(z)) if np.isfinite(z) else float("nan")
        rows.append((*c, est, se, z, p))
    return pd.DataFrame(rows, columns=[*cols, "est", "se", "z", "p"])


def sample_ploidy(prefix: str | os.PathLike[str], populations: Sequence[str],
                  pseudohaploid: str = "detect") -> dict[str, int]:
    """Ploidy of each sample of ``populations``: pseudo-haploid (1) or diploid (2).

    ``pseudohaploid="detect"``: pseudo-haploid when the sample has no
    heterozygous genotype at any SNP; ``"admixtools2"``: none among the first
    1000 SNPs, as admixtools2 assigns it (ADMIXTOOLS2-0007); ``"diploid"``.
    """
    return dict(_core.sample_ploidy(os.fspath(prefix), list(populations), pseudohaploid))


def _f3_rows(data: Data, combs: list[tuple], *, outgroupmode: bool, apply_corr: bool,
             block_size: float, auto_only: bool, pseudohaploid: str = "detect") -> pd.DataFrame:
    cols = ["pop1", "pop2", "pop3"]
    if _is_prefix(data):
        pops, triples = _index(combs)
        out = _core.f3_from_geno(os.fspath(data), pops, triples, outgroupmode, apply_corr, block_size,
                                 auto_only, pseudohaploid)
        return _stats_frame(combs, cols, out)
    return _from_f2(combs, cols, lambda c: data.f3(*c))


def outgroup_f3_matrix(
    data: Data,
    outgroup: str,
    populations: Sequence[str],
    *,
    value: str = "est",
    apply_corr: bool = True,
    block_size: float = 0.05,
    auto_only: bool = True,
    pseudohaploid: str = "detect",
) -> pd.DataFrame:
    """Outgroup f3(outgroup; A, B) for every pair of ``populations``.

    From genotypes this is ``qp3pop(outgroupmode=True)``; from f2 blocks, f3
    from the blocks. The diagonal is f3(outgroup; A, A), the drift from the
    outgroup to A.

    Parameters
    ----------
    data
        A genotype prefix, or f2 blocks.
    outgroup
        The f3 target, an outgroup to all of ``populations``.
    populations
        Rows and columns, in this order (see ``cluster_order``).
    value
        "est", "se", "z" or "p".

    Returns
    -------
    pandas.DataFrame
        Symmetric matrix indexed by ``populations``.
    """
    pops = list(populations)
    pairs = [(a, b) for i, a in enumerate(pops) for b in pops[i:]]
    rows = _f3_rows(data, [(outgroup, a, b) for a, b in pairs], outgroupmode=True, apply_corr=apply_corr,
                    block_size=block_size, auto_only=auto_only, pseudohaploid=pseudohaploid)
    m = pd.DataFrame(np.nan, index=pops, columns=pops)
    for (a, b), v in zip(pairs, rows[value]):  # noqa: B905 (strict= needs Python 3.10)
        m.loc[a, b] = m.loc[b, a] = v
    return m


def dstat_matrix(
    data: Data,
    populations: Sequence[str],
    *,
    pop1: str | None = None,
    pop2: str | None = None,
    pop3: str | None = None,
    pop4: str | None = None,
    value: str = "est",
    f4mode: bool = False,
    block_size: float = 0.05,
    auto_only: bool = True,
    pseudohaploid: str = "detect",
) -> pd.DataFrame:
    """D- (or f4-) statistics for every pair X, Y of ``populations``.

    Give one fixed pair: with ``pop3`` and ``pop4``, entry (X, Y) is
    D(X, Y; pop3, pop4); with ``pop1`` and ``pop2``, it is D(pop1, pop2; X, Y).
    Swapping X and Y negates the statistic (and z), so the matrix is
    antisymmetric, with NaN on the diagonal.

    Parameters
    ----------
    data
        A genotype prefix, or f2 blocks (f4 only).
    populations
        Rows and columns, in this order.
    pop1, pop2, pop3, pop4
        The fixed pair.
    value
        "est", "se", "z" or "p".
    f4mode
        True for f4 instead of D.

    Returns
    -------
    pandas.DataFrame
    """
    if (pop1 is None) == (pop3 is None) or (pop1 is None) != (pop2 is None) or (pop3 is None) != (pop4 is None):
        raise ValueError("give pop1 and pop2, or pop3 and pop4")
    pops = list(populations)
    pairs = [(a, b) for i, a in enumerate(pops) for b in pops[i + 1:]]
    combs = [(a, b, pop3, pop4) if pop3 is not None else (pop1, pop2, a, b) for a, b in pairs]
    if _is_prefix(data):
        names, quads = _index(combs)
        out = _core.f4_from_geno(os.fspath(data), names, quads, f4mode, block_size, auto_only, pseudohaploid)
        rows = _stats_frame(combs, ["pop1", "pop2", "pop3", "pop4"], out)
    else:
        if not f4mode:
            raise ValueError("D-statistics (f4mode=False) need genotypes, not f2 blocks")
        rows = _from_f2(combs, ["pop1", "pop2", "pop3", "pop4"], lambda c: data.f4(*c))
    sign = -1.0 if value in ("est", "z") else 1.0
    m = pd.DataFrame(np.nan, index=pops, columns=pops)
    for (a, b), v in zip(pairs, rows[value]):  # noqa: B905 (strict= needs Python 3.10)
        m.loc[a, b] = v
        m.loc[b, a] = sign * v
    return m


def cluster_order(matrix: pd.DataFrame) -> list[str]:
    """Row order that puts populations with similar rows next to each other.

    Average-linkage (UPGMA) clustering of the rows by Euclidean distance,
    missing values taken as 0; the leaves are read off the tree with the
    earlier-listed cluster first at each join, so the order is deterministic.

    Parameters
    ----------
    matrix
        A matrix with population rows (such as ``outgroup_f3_matrix``).

    Returns
    -------
    list of str
        The row labels in clustered order.
    """
    x = np.nan_to_num(matrix.to_numpy(dtype=float), nan=0.0)
    labels = list(matrix.index)
    d = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(axis=2))
    clusters: list[list[int]] = [[i] for i in range(len(labels))]
    while len(clusters) > 1:
        best = None
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                dist = d[np.ix_(clusters[i], clusters[j])].mean()
                if best is None or dist < best[0]:
                    best = (dist, i, j)
        _, i, j = best
        clusters[i] = clusters[i] + clusters[j]
        del clusters[j]
    return [labels[k] for k in clusters[0]] if clusters else []
