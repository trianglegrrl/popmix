"""Site-frequency spectra: admixtools2's joint SFS and dadi-style projected spectra."""

from __future__ import annotations

import itertools
import os
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from popmix import _core


@dataclass
class Sfs:
    """A site-frequency spectrum, as a dadi ``Spectrum``.

    Attributes
    ----------
    pops
        Populations, one axis each.
    data
        SNP counts (expected counts after projection), indexed by the count of
        the counted allele in each population.
    mask
        dadi's mask: the corners (allele absent or fixed everywhere) and,
        folded, the entries folded away.
    folded
        Whether the spectrum is folded (minor allele counts).
    not_integer
        SNPs left out because a count is half an allele (a heterozygous call
        in a pseudo-haploid sample).
    """

    pops: list[str]
    data: np.ndarray
    mask: np.ndarray
    folded: bool = False
    not_integer: int = 0

    @property
    def sample_sizes(self) -> list[int]:
        return [d - 1 for d in self.data.shape]

    def _summary(self) -> tuple[float, float, float, float]:
        if self.data.ndim != 1:
            raise ValueError("only defined for a one-population spectrum")
        if self.folded:
            raise ValueError("pi and Tajima's D need an unfolded spectrum")
        return _core.sfs_summary(self.data.tolist(), self.mask.tolist())

    def S(self) -> float:  # noqa: N802 (dadi's name)
        """Segregating sites.

        Returns
        -------
        float
            Number of segregating sites in the spectrum.
        """
        return self._summary()[0]

    def pi(self) -> float:
        """Expected pairwise differences (per SNP set), as dadi ``pi()``.

        Returns
        -------
        float
            Mean pairwise differences over the SNP set.
        """
        return self._summary()[1]

    def watterson_theta(self) -> float:
        """Watterson's theta, S over the harmonic number.

        Returns
        -------
        float
            ``S()`` divided by the sum of 1/i for i from 1 to n - 1 (n
            alleles).
        """
        return self._summary()[2]

    def tajima_d(self) -> float:
        """Tajima's D.

        Returns
        -------
        float
            Tajima's D, from ``pi()`` and ``watterson_theta()``.
        """
        return self._summary()[3]


def sfs(
    prefix: str | os.PathLike[str],
    populations: Sequence[str],
    projection: Sequence[int] | None = None,
    *,
    folded: bool = False,
    auto_only: bool = True,
    pseudohaploid: str = "detect",
) -> Sfs:
    """Site-frequency spectrum projected to ``projection`` alleles per population, as dadi.

    Each SNP is projected down to the given number of alleles in each
    population by sampling without replacement (hypergeometric), as dadi
    ``Spectrum.from_data_dict``: SNPs with missing data count in full as long
    as each population has at least the projected number of alleles
    observed, and are left out otherwise. Pseudo-haploid samples count one
    allele. Without ``projection``, the full sample sizes (SNPs with any
    missing data are left out).

    The counted allele is the one the genotype file counts (EIGENSTRAT's
    first allele, PLINK's A1), so the unfolded spectrum is polarized only if
    that allele is the derived one; otherwise use ``folded=True``.

    Parameters
    ----------
    prefix
        Genotype prefix.
    populations
        One axis each.
    projection
        Alleles per population to project to.
    folded
        Fold as dadi ``Spectrum.fold()``.
    auto_only
        Use only autosomes (the default: on X and Y, males carry one allele,
        which this does not take into account).
    pseudohaploid
        Which samples count one allele per SNP: "detect" (those with no
        heterozygous call at any SNP), "admixtools2" (no heterozygous call
        among the first 1000 SNPs, which can misclassify a diploid sample with
        missing data there) or "diploid" (none).

    Returns
    -------
    Sfs
    """
    pops = list(populations)
    r = _core.projected_sfs(os.fspath(prefix), pops, list(projection or []), folded, auto_only, pseudohaploid)
    return Sfs(pops=pops, data=np.asarray(r["data"]), mask=np.asarray(r["mask"]), folded=folded,
               not_integer=r["not_integer"])


def joint_sfs(
    prefix: str | os.PathLike[str],
    populations: Sequence[str],
    *,
    compat: str | None = None,
    auto_only: bool = False,
) -> pd.DataFrame:
    """Joint site-frequency spectrum as admixtools2 ``joint_sfs()``.

    For every combination of allele counts, from 0 to each population's
    largest sample size, the number of SNPs with those counts (the first
    population varies slowest). A SNP with missing data enters with its count
    among the observed alleles; a SNP where a population has no data is left
    out. For a spectrum that handles missing data properly, use ``sfs`` with a
    projection.

    Parameters
    ----------
    prefix
        Genotype prefix.
    populations
        Populations, one column each.
    compat
        None for popmix's correction, "admixtools2-2.0.10" to reproduce
        admixtools2 exactly (ADMIXTOOLS2-0003: it drops SNPs whose allele
        count is not exactly an integer in floating point).
    auto_only
        Use only autosomes. admixtools2 uses every SNP, counting every
        sample on X and Y with its autosomal ploidy.

    Returns
    -------
    pandas.DataFrame
        One column per population and ``n``. ``df.attrs`` holds n_snps,
        in_sfs, missing and not_integer.
    """
    pops = list(populations)
    r = _core.joint_sfs(os.fspath(prefix), pops, compat, auto_only)
    data = np.asarray(r["data"])
    grid = list(itertools.product(*(range(d) for d in data.shape)))
    df = pd.DataFrame(grid, columns=pops)
    df["n"] = data.ravel().astype(np.int64)
    df.attrs.update({k: r[k] for k in ("n_snps", "in_sfs", "missing", "not_integer")})
    return df
