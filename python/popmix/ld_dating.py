"""Admixture dating by the decay of admixture LD: DATES."""

from __future__ import annotations

import os
from collections.abc import Sequence

from popmix import _core
from popmix.results import DatesResult, dates_result


def dates(prefix: str | os.PathLike[str], admixed: str, sources: Sequence[str], *,
          generation_time: float = 29.0, generation_time_se: float = 0.0, sample_age: float = 0.0,
          sample_age_se: float = 0.0, min_distance_cm: float = 0.45,
          max_distance_cm: float = 100.0, bin_width_cm: float = 0.1, jackknife: bool = True,
          qbin: int = 0, compat: str | None = None) -> DatesResult:
    """Date the admixture of ``admixed`` from its two ``sources`` with DATES
    (Narasimhan et al. 2019): the decay of ancestry covariance with genetic
    distance.

    Parameters
    ----------
    prefix
        The genotype dataset; only the three populations' samples are read,
        their ploidy detected as for every analysis.
    admixed
        The admixed population (DATES's test population).
    sources
        Its two source populations.
    generation_time
        Years per generation, for ``years_bp``: an assumption (often 25 to
        30 years), so give its uncertainty as ``generation_time_se``.
    generation_time_se
        Standard error of the generation time, in years.
    sample_age
        Age of the admixed samples in years BP. DATES counts generations
        before the samples lived, so without it (0, present-day samples)
        ``years_bp`` is years before them, not before present.
    sample_age_se
        Standard error of the sample age, in years.
    min_distance_cm
        Fit bins from this distance (DATES ``lovalfit``).
    max_distance_cm
        Largest distance (DATES ``maxdis`` x 100).
    bin_width_cm
        Distance bin width (DATES ``binsize`` x 100).
    jackknife
        Leave-one-chromosome-out jackknife for the SE.
    qbin
        DATES ``qbin``: 0 visits every SNP pair; N > 0 uses an FFT mesh of N
        cells per bin, faster and slightly approximate.
    compat
        ``"dates-v753"`` reproduces DATES v753's known quirks, for output
        identical to it (DEVIATIONS.md); None corrects them.

    Returns
    -------
    DatesResult
        Generations since admixture with SE, the LD bins and the fit.
    """
    return dates_result(_core.dates(
        os.fspath(prefix), admixed, list(sources), generation_time=generation_time,
        generation_time_se=generation_time_se, sample_age=sample_age, sample_age_se=sample_age_se,
        min_distance_cm=min_distance_cm, max_distance_cm=max_distance_cm,
        bin_width_cm=bin_width_cm, jackknife=jackknife, qbin=qbin, compat=compat))
