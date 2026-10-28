"""Each sample's coverage: SNPs called, heterozygosity and ploidy (TASKS.md H10).

Ancient samples are often of low coverage: most SNPs uncalled, and one
allele read at the rest (pseudo-haploid, coded 0 or 2). The table says
which samples to project rather than use to build PCA axes, which share
too few SNPs for KING-robust (use READ), and on how many calls each
ploidy rule decides.
"""

from __future__ import annotations

import os
from collections.abc import Sequence

import pandas as pd

from popmix import _core


def sample_coverage(prefix: str | os.PathLike[str], populations: Sequence[str] | None = None, *,
                    auto_only: bool = True) -> pd.DataFrame:
    """Each sample's calls and ploidy.

    Parameters
    ----------
    prefix
        Genotype prefix.
    populations
        Populations whose samples are listed; all (but Ignore) if None.
    auto_only
        Count only autosomal SNPs.

    Returns
    -------
    pandas.DataFrame
        One row per sample, in file order: ``id``, ``pop``, ``snps_called``,
        ``called_fraction`` (of the SNPs read), ``heterozygosity`` (the share
        of calls that are heterozygous; 0 for a pseudo-haploid sample),
        ``ploidy`` (1 pseudo-haploid or 2, as popmix finds it: no
        heterozygous call at any SNP) and ``ploidy_admixtools2`` (none among
        the first 1,000 SNPs, as admixtools2).
    """
    prefix = os.fspath(prefix)
    pops = list(populations) if populations is not None else \
        [p for p, _ in _core.list_populations(prefix) if p != "Ignore"]
    rows, n_snps = _core.sample_calls(prefix, pops, auto_only=auto_only)
    df = pd.DataFrame(rows, columns=["id", "pop", "snps_called", "heterozygous"])
    df["called_fraction"] = df["snps_called"] / n_snps if n_snps else float("nan")
    df["heterozygosity"] = df["heterozygous"] / df["snps_called"].where(df["snps_called"] > 0)
    for column, rule in (("ploidy", "detect"), ("ploidy_admixtools2", "admixtools2")):
        ploidy = dict(_core.sample_ploidy(prefix, pops, rule))
        df[column] = df["id"].map(ploidy).astype("Int64")
    return df[["id", "pop", "snps_called", "called_fraction", "heterozygosity", "ploidy", "ploidy_admixtools2"]]
