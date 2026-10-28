"""Kinship between samples."""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import TYPE_CHECKING

from popmix import _core
from popmix._names import names

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd


def king(prefix: str | os.PathLike[str], populations: Sequence[str] | None = None) -> pd.DataFrame:
    """KING-robust kinship for every pair of samples, as PLINK 2 ``--make-king-table``.

    The between-family estimator of Manichaikul et al. (2010), from the SNPs
    called in both samples: 1/2 - (4 IBS0 + HET1_HOM2 + HET2_HOM1) /
    (4 min(het1, het2)). About 0.5 for duplicates or monozygotic twins, 0.25
    first-degree, 0.125 second-degree, 0 unrelated. It needs diploid calls: a
    sample with no heterozygous call (such as a pseudo-haploid ancient
    sample) gives minus infinity, as in PLINK 2, with a ``PopmixWarning``.
    Genotyping error at low coverage and contamination (false heterozygotes)
    bias it.

    Parameters
    ----------
    prefix
        Genotype prefix (autosomes are used).
    populations
        Populations whose samples are compared; all if None.

    Returns
    -------
    pandas.DataFrame
        One row per pair, in PLINK 2's order ((2nd, 1st), (3rd, 1st), (3rd,
        2nd), ...): id1, id2, nsnp, hethet, ibs0, het1_hom2, het2_hom1 (counts)
        and kinship.
    """
    import pandas as pd

    return pd.DataFrame(king_columns(prefix, populations))


def king_columns(
    prefix: str | os.PathLike[str], populations: Sequence[str] | None = None
) -> dict[str, list[str] | np.ndarray]:
    """The columns of :func:`king`, as a dict of lists and arrays (no pandas)."""
    if populations is None:
        populations = [p for p, _ in _core.list_populations(os.fspath(prefix)) if p != "Ignore"]
    return _core.king(os.fspath(prefix), names(populations))


def king_tsv(prefix: str | os.PathLike[str], populations: Sequence[str] | None = None) -> str:
    """The table of :func:`king` as TSV text, built in Rust (the kinship
    command's fast path: no numpy, no loop over the pairs in Python)."""
    if populations is None:
        populations = [p for p, _ in _core.list_populations(os.fspath(prefix)) if p != "Ignore"]
    return _core.king_tsv(os.fspath(prefix), names(populations))


def read_kinship(
    prefix: str | os.PathLike[str],
    populations: Sequence[str] | None = None,
    *,
    norm: str = "median",
    norm_value: float | None = None,
    window_size: int = 5_000_000,
    window_estimate: bool = False,
    thresholds: str = "read",
    compat: str | None = None,
    auto_only: bool = False,
) -> pd.DataFrame:
    """READ-style kinship for pseudo-haploid (ancient) samples, as READv2.

    For each pair, P0 is the proportion of SNPs called in both where they
    differ (a heterozygote counts as half), normalized by the P0 of unrelated
    pairs: by default the median over all pairs, which assumes most pairs are
    unrelated (use ``norm="max"`` for a trio, or ``norm_value`` from another
    sample of the population). Normalized P0 is about 1 for unrelated pairs,
    0.75 first degree, 0.875 second, 0.9375 third and 0.5 for the same
    individual; the kinship coefficient is 1 minus it. The class thresholds
    are fixed, so a norm from fewer than 10 pairs warns; background
    relatedness or inbreeding shifts every class, and contamination raises P0
    (relatives look more distant).

    Parameters
    ----------
    prefix
        Genotype prefix.
    populations
        Populations whose samples are compared (as one group, so they share
        the normalization); all if None.
    norm
        "median", "mean" or "max" of P0 over the pairs.
    norm_value
        Normalize by this P0 instead.
    window_size
        Jackknife window in base pairs (READv2 ``--window_size``).
    window_estimate
        P0 as the mean over windows and its SE across them (READv2
        ``--window_est``, READv1's method).
    thresholds
        "read" (READv2's) or "2pow" (READv2 ``--2pow``).
    compat
        None for popmix's corrections, "readv2-2.01" to reproduce READv2
        exactly (READV2-0001 to 0003 in DEVIATIONS.md).
    auto_only
        Use only autosomes (READv2 uses every SNP in its input).

    Returns
    -------
    pandas.DataFrame
        One row per pair (pairs with no SNP called in both are left out):
        id1, id2, relationship (READv2's labels), z_upper, z_lower,
        normalized_p0, normalized_p0_se, p0, p0_se, first_degree_type
        ("Parent-offspring", "Siblings", "NotApplicable" or ""),
        first_degree_windows_outside, overlap_snps, effective_snps, kinship,
        kinship_ci_low, kinship_ci_high. ``df.attrs["norm_value"]`` holds the
        normalization value.
    """
    if populations is None:
        populations = [p for p, _ in _core.list_populations(os.fspath(prefix)) if p != "Ignore"]
    out = _core.read_kinship(os.fspath(prefix), names(populations), norm, norm_value, window_size,
                             window_estimate, thresholds, compat, auto_only)
    norm_used = out.pop("norm_value")
    import pandas as pd

    df = pd.DataFrame(out)
    df.attrs["norm_value"] = norm_used
    return df
