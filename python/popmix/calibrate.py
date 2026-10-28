"""Calibration by simulation: how often a design rejects (TASKS.md H7).

Run one qpAdm model, or one f4 statistic, on each of many simulated
datasets and count the rejections at ``alpha``. When the simulations make
the model true (a qpAdm model with the right sources; an f4 whose four
populations form a tree), the rejection rate is the design's false-positive
rate, which should be near ``alpha``. When they make it false, it is its
power. The rate comes with a Wilson 95% interval, so a handful of
replicates is not over-read.

popmix reads the simulated genotypes; it does not simulate. Make them with
msprime (or any simulator) in your own script, write each replicate as a
VCF or an EIGENSTRAT/PLINK dataset, and pass the files here. msprime is
GPL-3 and popmix (MIT) never imports it.
"""

from __future__ import annotations

import math
import os
from collections.abc import Sequence
from typing import Any

import pandas as pd

from popmix import _core
from popmix.fits import qpadm


def wilson_interval(rejected: int, n: int, confidence: float = 0.95) -> tuple[float, float]:
    """Wilson score interval of a proportion ``rejected / n``.

    Parameters
    ----------
    rejected, n
        Successes and trials.
    confidence
        Coverage of the interval (0.95: z = 1.96).

    Returns
    -------
    tuple of float
        Lower and upper bounds, within [0, 1]; (0, 1) for ``n == 0``.
    """
    if n == 0:
        return 0.0, 1.0
    from statistics import NormalDist

    z = NormalDist().inv_cdf(0.5 + confidence / 2)
    p = rejected / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, centre - half), min(1.0, centre + half)


def _summarize(df: pd.DataFrame, alpha: float) -> pd.DataFrame:
    n, k = len(df), int(df["rejected"].sum())
    lo, hi = wilson_interval(k, n)
    df.attrs.update({"alpha": alpha, "n": n, "rejected": k, "rate": k / n if n else float("nan"),
                     "ci_low": lo, "ci_high": hi})
    return df


def _f2(replicate: Any, pops: list[str], f2_options: dict) -> _core.F2Blocks:
    if isinstance(replicate, (str, os.PathLike)):
        return _core.extract_f2(os.fspath(replicate), pops, **f2_options)
    return replicate


def calibrate_qpadm(replicates: Sequence[Any], left: Sequence[str], right: Sequence[str], target: str, *,
                    alpha: float = 0.05, constrained: bool = False, **f2_options: Any) -> pd.DataFrame:
    """One qpAdm model on each replicate, and how often it is rejected.

    Parameters
    ----------
    replicates
        Genotype prefixes (or VCF paths), or f2 blocks, one per replicate.
    left, right, target
        The model, as for ``qpadm``.
    alpha
        A replicate rejects the model when p is at or below it.
    constrained
        As for ``qpadm``.
    **f2_options
        Passed to ``extract_f2`` for each prefix (``block_size``,
        ``pseudohaploid``, ``inds`` and ``pops`` to label VCF samples, ...).

    Returns
    -------
    pandas.DataFrame
        One row per replicate: ``replicate``, ``chisq``, ``dof``, ``p``,
        ``feasible``, ``weights`` (comma-separated) and ``rejected``. Its
        ``attrs`` hold ``n``, ``rejected``, ``rate`` and the Wilson 95%
        interval ``ci_low``, ``ci_high``.
    """
    pops = list(dict.fromkeys([target, *left, *right]))
    rows = []
    for k, rep in enumerate(replicates):
        r = qpadm(_f2(rep, pops, f2_options), left, right, target, constrained=constrained)
        rows.append({"replicate": os.fspath(rep) if isinstance(rep, (str, os.PathLike)) else k,
                     "chisq": r.chisq, "dof": r.dof, "p": r.p, "feasible": r.feasible,
                     "weights": ",".join(f"{w:.6g}" for w in r.weights), "rejected": bool(r.p <= alpha)})
    return _summarize(pd.DataFrame(rows), alpha)


def calibrate_f4(replicates: Sequence[Any], pop1: str, pop2: str, pop3: str, pop4: str, *,
                 alpha: float = 0.05, **f2_options: Any) -> pd.DataFrame:
    """f4(pop1, pop2; pop3, pop4) on each replicate, and how often it is
    significantly non-zero.

    Parameters
    ----------
    replicates
        Genotype prefixes (or VCF paths), or f2 blocks, one per replicate.
    pop1, pop2, pop3, pop4
        The statistic.
    alpha
        A replicate rejects f4 = 0 when its two-sided p is at or below it.
    **f2_options
        Passed to ``extract_f2`` for each prefix.

    Returns
    -------
    pandas.DataFrame
        One row per replicate: ``replicate``, ``est``, ``se``, ``z``, ``p``
        and ``rejected``, with ``attrs`` as ``calibrate_qpadm``'s.
    """
    pops = list(dict.fromkeys([pop1, pop2, pop3, pop4]))
    rows = []
    for k, rep in enumerate(replicates):
        r = _f2(rep, pops, f2_options).f4(pop1, pop2, pop3, pop4)
        p = r["p"]
        rows.append({"replicate": os.fspath(rep) if isinstance(rep, (str, os.PathLike)) else k,
                     "est": r["est"], "se": r["se"], "z": r["z"], "p": p,
                     "rejected": bool(p <= alpha)})
    return _summarize(pd.DataFrame(rows), alpha)
