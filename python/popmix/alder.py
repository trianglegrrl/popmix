"""ALDER: weighted LD curves for dating and testing admixture (Loh et al. 2013)."""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from popmix import _core


@dataclass
class AlderCurve:
    """One weighted LD curve and its fits.

    Attributes
    ----------
    refs
        The references whose weights were used.
    fit_start
        Distance (Morgans) from which the curve is fit.
    curve
        d_cm, wld, count per distance bin; the last row (d_cm = inf) is the
        affine term from SNP pairs on different chromosomes, if computed.
    fits
        One row per fit start (fit_start - 0.2 cM to + 0.2 cM): mindis_cm,
        decay (generations), amp_exp, amp_tot, amp_aff, each with its
        jackknife SE, and p (from the smaller z-score of decay and amp_exp).
    fit_index
        The row of ``fits`` fit from ``fit_start``.
    """

    refs: list[str]
    fit_start: float
    curve: pd.DataFrame
    fits: pd.DataFrame
    fit_index: int
    jackknife_curves: list[pd.DataFrame] = field(default_factory=list)

    @property
    def fit(self) -> pd.Series:
        return self.fits.iloc[self.fit_index]


@dataclass
class AlderResult:
    """Result of ``alder``.

    Attributes
    ----------
    test_pop, refs, chromosomes
        The test population, the references, and the chromosomes (jackknife
        replicates).
    ld_corr_stops
        Per reference, the distance (Morgans) to which the test population's
        LD is correlated with its LD (infinite past 2 cM); None if ``mindis``
        was given.
    curve
        The weighted LD curve with the references' weights (two references:
        their allele frequency difference).
    oneref_curves
        With two references, the one-reference curve with each.
    test
        With two references, the admixture test: pretest (both one-reference
        curves significant), tworef_significant, consistent (decay rates
        within 25%), succeeds, z, p and decay_diffs.
    mixture_fraction
        With one reference, the lower bound on the mixture fraction (value,
        SE).
    """

    test_pop: str
    refs: list[str]
    chromosomes: list[str]
    ld_corr_stops: dict
    curve: AlderCurve
    oneref_curves: dict = field(default_factory=dict)
    test: dict | None = None
    mixture_fraction: tuple | None = None


def _frame(c: dict) -> pd.DataFrame:
    return pd.DataFrame({"d_cm": 100 * np.asarray(c["d"]), "wld": c["wld"], "count": c["count"]}).iloc[1:]


def _curve(run: dict, refs: list[str]) -> AlderCurve:
    rows = []
    for f in run["fits"]:
        row = {"mindis_cm": 100 * f["mindis"]}
        for var in ("decay", "amp_exp", "amp_tot", "amp_aff"):
            row[var] = f[var]["value"]
            row[f"{var}_se"] = f[var]["se"]
        row["p"] = f["p"]
        row["inter_chrom_affine"] = f["inter_chrom_affine"]
        row["jackknife_failures"] = f["failures"]
        rows.append(row)
    return AlderCurve(refs=[refs[i] for i in run["refs"]], fit_start=run["fit_start"], curve=_frame(run["curve"]),
                      fits=pd.DataFrame(rows), fit_index=run["fit_index"],
                      jackknife_curves=[_frame(c) for c in run["jackknife_curves"]])


def alder(
    prefix: str | os.PathLike[str],
    test: str,
    refs: Sequence[str],
    *,
    mindis: float | None = None,
    binsize: float = 0.0005,
    maxdis: float = 0.5,
    mincount: int = 4,
    compat: str | None = None,
) -> AlderResult:
    """ALDER weighted LD for ``test`` with one or two reference populations.

    With two references A and B, the test population's LD is weighted by the
    A - B allele frequency difference, fit as A exp(-n d) + C (n the
    generations since admixture), and checked with ALDER's test: the
    two-reference curve and both one-reference curves (each reference with
    the test population) must be significant. With one reference, the curve
    and a lower bound on the mixture fraction.

    The curve is fit from where the test population's LD stops being
    correlated with each reference's (at least 0.5 cM) unless ``mindis`` is
    given; standard errors are from a jackknife over chromosomes. Autosomes
    1-22 are used, sorted by genetic position; the .snp file needs a real
    genetic map.

    Parameters
    ----------
    prefix
        Genotype prefix.
    test
        Test (possibly admixed) population.
    refs
        One or two reference populations.
    mindis
        Fit start (Morgans), skipping the LD correlation step.
    binsize, maxdis
        Distance bin and largest distance (Morgans).
    mincount
        Test samples that must be called at a SNP (two-reference curve).
    compat
        None for popmix's correction, "alder-1.03" to reproduce ALDER exactly
        (ALDER-0001 in DEVIATIONS.md).

    Returns
    -------
    AlderResult
    """
    refs = list(refs)
    r = _core.alder(os.fspath(prefix), test, refs, mindis, binsize, maxdis, mincount, compat)
    return AlderResult(
        test_pop=test,
        refs=refs,
        chromosomes=list(r["chromosomes"]),
        ld_corr_stops=dict(zip(refs, r["ld_corr_stops"])),  # noqa: B905
        curve=_curve(r["run"], refs),
        oneref_curves={refs[i]: _curve(run, refs) for i, run in enumerate(r["oneref_runs"])},
        test=r["test"],
        mixture_fraction=None if r["mixture_fraction"] is None else tuple(r["mixture_fraction"]),
    )
