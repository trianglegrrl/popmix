"""ROLLOFF: admixture dates from the decay of weighted LD correlation (ADMIXTOOLS)."""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from popmix import _core


@dataclass
class RolloffResult:
    """Result of ``rolloff``.

    Attributes
    ----------
    bins
        One row per distance bin, labeled by its upper edge as ROLLOFF does:
        d_cm, zscore (sum(w z) / sqrt(sum w^2)), regression, corr (the
        correlation of weights and LD, which is fit) and count (SNP pairs).
    decay, affine, amplitude
        The fit corr = C + A exp(-m d) from ``lval`` to ``hval`` cM: m in
        generations, C and A.
    decay_jackknife
        Weighted jackknife over chromosomes (mean, standard error).
    jackknife_bins
        The bins with each chromosome left out.
    """

    bins: pd.DataFrame
    decay: float
    affine: float
    amplitude: float
    decay_jackknife: tuple
    chromosomes: list[str]
    n_snps_used: int
    jackknife_bins: list[pd.DataFrame] = field(default_factory=list)


def _frame(b: dict) -> pd.DataFrame:
    return pd.DataFrame({"d_cm": 100 * np.asarray(b["d"]), "zscore": b["zscore"], "regression": b["regression"],
                         "corr": b["corr"], "count": b["count"]})


def rolloff(
    prefix: str | os.PathLike[str],
    admixed: str,
    refs: Sequence[str],
    *,
    binsize: float = 0.0005,
    maxdis: float = 0.05,
    runmode: int = 1,
    zdipcorr: bool = True,
    chithresh: float = 6.0,
    mincount: int = 5,
    minparentcount: int = 10,
    chromosomes: Sequence[str | int] | None = None,
    lval: float = 0.5,
    hval: float = 100.0,
) -> RolloffResult:
    """ROLLOFF (ADMIXTOOLS 7.0.2) for ``admixed`` with two reference populations.

    For SNP pairs closer than ``maxdis``, the admixed population's LD (with
    ``zdipcorr``, its Fisher-transformed genotype correlation) is correlated
    with the product of the references' allele frequency differences, bin
    by bin, leaving out pairs in LD in either reference (chi-square at least
    ``chithresh``). The date is the decay rate of C + A exp(-m d) fit to
    the correlations (least squares, as rexpfit.r's nls), with a weighted
    jackknife over chromosomes (weights: SNPs per chromosome in the .snp
    file, as expfit.sh). DATES (``popmix.dates``) and ALDER (``alder``) are
    its successors.

    Parameters
    ----------
    prefix
        Genotype prefix with a genetic map.
    admixed
        The admixed population.
    refs
        The two reference populations.
    binsize, maxdis
        Distance bin and largest distance (Morgans).
    runmode
        1 (weights pA - pB, recommended) or 0 (divided by sqrt of the mean
        frequency).
    zdipcorr, chithresh, mincount, minparentcount
        ROLLOFF's ``zdipcorrmode``, ``chithresh``, ``mincount`` and
        ``minparentcount``.
    chromosomes
        Use only these chromosomes.
    lval, hval
        Fit range (cM).

    Returns
    -------
    RolloffResult
    """
    chroms = None if chromosomes is None else [str(c) for c in chromosomes]
    r = _core.rolloff(os.fspath(prefix), admixed, list(refs), binsize, maxdis, runmode, zdipcorr, chithresh,
                      mincount, minparentcount, chroms, lval, hval)
    c, a, m = r["fit"]
    return RolloffResult(bins=_frame(r["bins"]), decay=m, affine=c, amplitude=a, decay_jackknife=tuple(r["jackknife"]),
                         chromosomes=list(r["chromosomes"]), n_snps_used=r["n_snps_used"],
                         jackknife_bins=[_frame(b) for b in r["jackknife_bins"]])
