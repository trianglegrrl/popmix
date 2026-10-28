"""PCA as EIGENSOFT's smartpca, with projection of other samples."""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from popmix import _core
from popmix.pca_holdout import measured_shrinkage, warn_shrinkage


@dataclass
class PcaResult:
    """Result of ``smartpca``.

    Attributes
    ----------
    eigenvalues
        Every eigenvalue of the scaled genotype matrix, decreasing.
    coords
        One row per sample with coordinates (PCA samples and projected
        samples, in .ind order): id, pop, role ("pca", "project", or with
        ``holdout`` "holdout") and PC1, PC2, ...
    tracy_widom
        eigenvalue, twstat, p and neff (effective number of markers) for each
        eigenvalue; NaN where too few eigenvalues are left for a statistic.
    outliers
        Samples removed as outliers: id, pop, iteration (from 1), pc (from
        1) and score (standard deviations).
    no_data
        Samples left out for lack of data.
    shrinkage
        With ``holdout``, each component's measured shrinkage: the slope of
        the held-out samples' population means on the PCA samples' (1 for
        none, NaN where the populations do not separate); empty without.
    missing_fraction
        Each sample's share of the PCA's SNPs it misses, in ``coords`` order.
    """

    eigenvalues: np.ndarray
    coords: pd.DataFrame
    tracy_widom: pd.DataFrame
    outliers: pd.DataFrame
    no_data: list[str] = field(default_factory=list)
    n_pca_samples: int = 0
    n_snps_used: int = 0
    n_snps_killed: int = 0
    shrinkage: np.ndarray = field(default_factory=lambda: np.zeros(0))
    missing_fraction: np.ndarray = field(default_factory=lambda: np.zeros(0))


def smartpca(
    prefix: str | os.PathLike[str],
    pca_pops: Sequence[str] | None = None,
    project_pops: Sequence[str] = (),
    *,
    n_components: int = 10,
    outlier_iterations: int = 5,
    outlier_components: int = 10,
    outlier_sigma: float = 6.0,
    lsqproject: bool = False,
    min_tw_eigs: int = 10,
    compat: str | None = None,
    holdout: bool = False,
) -> PcaResult:
    """PCA computed on ``pca_pops``, with ``project_pops`` projected onto it.

    The same computation as EIGENSOFT 7.2.1 smartpca with ``poplistname``
    (smartpca's defaults are the defaults here): autosomal SNPs polymorphic
    in the PCA samples, each centered and scaled by 1/sqrt(p(1-p)); the
    sample matrix scaled to trace n - 1; outliers more than
    ``outlier_sigma`` standard deviations out on one of the first
    ``outlier_components`` eigenvectors removed, up to
    ``outlier_iterations`` times; Tracy-Widom statistics for the
    eigenvalues. Projected samples are placed by the dot product of their
    normalized genotypes with the SNP loadings, or with ``lsqproject`` by
    least squares on their observed genotypes (recommended for ancient
    samples with missing data: by dot product a missing genotype counts as
    0, which pulls a sample toward the origin, and a projected sample
    missing more than 10% of the SNPs warns).

    Projected samples are also shrunk toward the origin relative to the
    PCA samples, whose own noise helped make the axes: two halves of one
    population, one in the PCA and one projected, do not coincide (on
    six 1000 Genomes populations, by 2% on PC1 to 72% on PC6). With
    ``holdout=True`` every second sample of each PCA population, in .ind
    order, is projected instead (role "holdout"): projected samples are
    then compared with reference samples placed the same way, and
    ``shrinkage`` measures each component's shrinkage, with a warning
    where it is below 0.9. The axes are made from half the samples.

    Parameters
    ----------
    prefix
        Genotype prefix.
    pca_pops
        Populations the PCA is computed on; all populations (except
        ``project_pops`` and "Ignore") if None.
    project_pops
        Populations projected onto the PCA.
    n_components, outlier_iterations, outlier_components, outlier_sigma
        smartpca ``numoutevec``, ``numoutlieriter``, ``numoutlierevec``,
        ``outliersigmathresh``.
    lsqproject
        smartpca ``lsqproject: YES``.
    min_tw_eigs
        smartpca ``nostatslim``.
    compat
        None for popmix's corrections, or "eigensoft-7.2.1" to reproduce
        smartpca exactly (EIGENSOFT-0001 in DEVIATIONS.md).
    holdout
        Hold out every second sample of each PCA population (see above);
        needs two PCA populations with at least 2 samples.

    Returns
    -------
    PcaResult
        Eigenvector signs are arbitrary, as in smartpca.
    """
    project = list(project_pops)
    if pca_pops is None:
        pca = [p for p, _ in _core.list_populations(os.fspath(prefix)) if p not in project and p != "Ignore"]
    else:
        pca = list(pca_pops)
    r = _core.smartpca(os.fspath(prefix), pca, project, n_components, outlier_iterations, outlier_components,
                       outlier_sigma, lsqproject, min_tw_eigs, compat, holdout=holdout)
    coords = np.asarray(r.coords)
    df = pd.DataFrame({"id": r.sample_ids, "pop": r.populations, "role": r.roles})
    for j in range(coords.shape[1] if coords.ndim == 2 else 0):
        df[f"PC{j + 1}"] = coords[:, j]
    tw = pd.DataFrame(r.tracy_widom, columns=["eigenvalue", "twstat", "p", "neff"]).astype(float)
    outliers = pd.DataFrame(r.outliers, columns=["id", "pop", "iteration", "pc", "score"])
    outliers["pc"] += 1
    shrinkage = measured_shrinkage(df) if holdout else np.zeros(0)
    if holdout:
        warn_shrinkage(shrinkage)
    return PcaResult(
        eigenvalues=np.asarray(r.eigenvalues),
        coords=df,
        tracy_widom=tw,
        outliers=outliers,
        no_data=list(r.no_data),
        n_pca_samples=r.n_pca_samples,
        n_snps_used=r.n_snps_used,
        n_snps_killed=r.n_snps_killed,
        shrinkage=shrinkage,
        missing_fraction=np.asarray(r.missing_fraction, dtype=float),
    )
