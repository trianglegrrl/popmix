"""Ancestry proportions as ADMIXTURE: unsupervised, supervised and projection."""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from popmix import _core
from popmix._names import names


@dataclass
class AdmixtureResult:
    """Result of ``admixture``.

    Attributes
    ----------
    Q
        One row per sample: id, pop and Q1..QK (ancestry proportions).
    P
        Allele frequencies, SNPs x K, of the allele the genotype files do not
        count (ADMIXTURE's .P).
    snps
        SNP IDs, the rows of P (SNPs with no genotype at all are left out).
    loglik
        Log-likelihood, as ADMIXTURE reports it.
    cv_error
        Cross-validation error, if requested.
    """

    Q: pd.DataFrame
    P: np.ndarray
    snps: list[str]
    loglik: float
    iterations: int
    converged: bool
    cv_error: float | None = None

    @property
    def k(self) -> int:
        return self.P.shape[1]

    def q_matrix(self) -> np.ndarray:
        """The ancestry proportions as a matrix.

        Returns
        -------
        numpy.ndarray
            Samples x K: columns Q1..QK of ``Q``, rows in the order of ``Q``.
        """
        return self.Q[[f"Q{j + 1}" for j in range(self.k)]].to_numpy()


def admixture(
    prefix: str | os.PathLike[str],
    k: int,
    populations: Sequence[str] | None = None,
    *,
    supervised: Sequence[str] | None = None,
    fixed_p: np.ndarray | None = None,
    cv: int | None = None,
    chromosomes: Sequence[str | int] | None = None,
    snps: Sequence[str] | None = None,
    seed: int = 1,
    tolerance: float = 1e-4,
    max_iterations: int = 10_000,
) -> AdmixtureResult:
    """Maximum-likelihood ancestry proportions with K ancestral populations.

    The model and likelihood are ADMIXTURE's (Alexander et al. 2009), with Q
    and P kept within [1e-5, 1 - 1e-5]; it is maximized by EM with SQUAREM
    acceleration from a random start. Genotypes are read as called, so a
    pseudo-haploid sample counts as homozygous, as in ADMIXTURE.

    Parameters
    ----------
    prefix
        Genotype prefix.
    k
        Number of ancestral populations.
    populations
        Populations whose samples are used; all if None.
    supervised
        Reference populations (K of them): their samples keep Q at their
        population (ADMIXTURE ``--supervised``). Clusters are numbered in the
        order the populations first appear, as ADMIXTURE numbers them.
    fixed_p
        Projection (ADMIXTURE ``-P``): P fixed to this (SNPs x K, the P of an
        earlier fit on the same SNPs); only Q is estimated.
    cv
        Folds of cross-validation (ADMIXTURE ``--cv``): the mean deviance of
        held-out genotypes, for choosing K.
    chromosomes
        Use only SNPs on these chromosomes.
    snps
        Use exactly these SNPs (IDs; in file order), for example the SNPs of
        the fit whose P is given as ``fixed_p``. Otherwise SNPs with no
        genotype in any sample are left out.
    seed
        Seed for the random start (and the folds).
    tolerance
        Stop when the log-likelihood improves by less than this (ADMIXTURE
        ``-C``).
    max_iterations
        Upper limit on accelerated iterations.

    Returns
    -------
    AdmixtureResult
    """
    if populations is None:
        populations = [p for p, _ in _core.list_populations(os.fspath(prefix)) if p != "Ignore"]
    chroms = None if chromosomes is None else [str(c) for c in chromosomes]
    p = None if fixed_p is None else np.ascontiguousarray(fixed_p, dtype=float)
    r = _core.admixture(os.fspath(prefix), k, names(populations), None if supervised is None else list(supervised),
                        p, cv, chroms, None if snps is None else list(snps), seed, tolerance, max_iterations)
    q = np.asarray(r["Q"])
    df = pd.DataFrame({"id": r["ids"], "pop": r["pops"]})
    for j in range(k):
        df[f"Q{j + 1}"] = q[:, j]
    return AdmixtureResult(Q=df, P=np.asarray(r["P"]), snps=list(r["snps"]), loglik=r["loglik"],
                           iterations=r["iterations"], converged=r["converged"], cv_error=r["cv_error"])


def admixture_loglik(
    prefix: str | os.PathLike[str],
    q: np.ndarray,
    p: np.ndarray,
    populations: Sequence[str] | None = None,
    *,
    chromosomes: Sequence[str | int] | None = None,
    snps: Sequence[str] | None = None,
) -> float:
    """ADMIXTURE's log-likelihood of given Q (samples x K) and P (SNPs x K).

    The samples and SNPs are those ``admixture`` would use with the same
    ``populations`` and ``chromosomes``.

    Parameters
    ----------
    prefix
        Genotype prefix.
    q
        Ancestry proportions, samples x K.
    p
        Allele frequencies, SNPs x K, as ``AdmixtureResult.P``.
    populations
        Populations whose samples are used; all if None.
    chromosomes
        Use only SNPs on these chromosomes.
    snps
        Use exactly these SNPs (IDs; in file order).

    Returns
    -------
    float
        The log-likelihood, as ADMIXTURE reports it.
    """
    if populations is None:
        populations = [pp for pp, _ in _core.list_populations(os.fspath(prefix)) if pp != "Ignore"]
    chroms = None if chromosomes is None else [str(c) for c in chromosomes]
    return _core.admixture_loglik(os.fspath(prefix), names(populations), np.ascontiguousarray(q, dtype=float),
                                  np.ascontiguousarray(p, dtype=float), chroms, None if snps is None else list(snps))
