"""f2 directories in admixtools2's format, and population subsets.

admixtools2 ``extract_f2(pref, outdir)`` writes, for every pair of
populations (pop1 <= pop2 in C order, each population with itself too),
``outdir/pop1/pop2_f2.rds``: a two-column matrix of each block's f2 and the
fraction of its SNPs with data for both, and ``outdir/block_lengths_f2.rds``
(``read_f2(counts = TRUE)`` multiplies the fraction by the block length). ``read_admixtools2_f2``
reads such a directory into f2 blocks and ``write_admixtools2_f2`` writes
one admixtools2 reads (``f2_from_precomp``, ``read_f2``).

admixtools2 also stores allele-frequency products (``*_ap.rds``) for
qpAdm on a directory, and FST; popmix's f2 blocks hold neither, so they are
not read or written. Give admixtools2's ``qpadm`` the f2 array
(``f2_from_precomp(dir)``) instead of the directory.
"""

from __future__ import annotations

import math
import os
from collections.abc import Sequence

from popmix import _core
from popmix.rds import matrix, read_rds, write_rds


def _c_order(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a.encode() <= b.encode() else (b, a)


def _vector(x) -> list:
    return x.data if hasattr(x, "data") else list(x)


def read_admixtools2_f2(f2_dir: str | os.PathLike, populations: Sequence[str] | None = None,
                        *, remove_na: bool = True) -> _core.PyF2Blocks:
    """Read an admixtools2 f2 directory (``extract_f2`` output).

    Parameters
    ----------
    f2_dir
        The directory.
    populations
        Populations to read (default: every subdirectory, in sorted order as
        admixtools2 lists them).
    remove_na
        Leave out blocks where some pair has no data, as admixtools2
        ``read_f2`` does.
    """
    f2_dir = os.fspath(f2_dir)
    pops = list(populations) if populations is not None else sorted(
        d for d in os.listdir(f2_dir) if os.path.isdir(os.path.join(f2_dir, d)))
    lengths_file = os.path.join(f2_dir, "block_lengths_f2.rds")
    if not os.path.exists(lengths_file):
        raise FileNotFoundError(f"{lengths_file} not found: not an admixtools2 f2 directory")
    lengths = [int(x) for x in _vector(read_rds(lengths_file))]
    nb = len(lengths)
    means = [[None] * len(pops) for _ in pops]
    counts = [[None] * len(pops) for _ in pops]
    for i, a in enumerate(pops):
        for j in range(i, len(pops)):
            b = pops[j]
            p1, p2 = _c_order(a, b)
            fl = os.path.join(f2_dir, p1, f"{p2}_f2.rds")
            if not os.path.exists(fl):
                fl = os.path.join(f2_dir, p2, f"{p1}_f2.rds")
            if not os.path.exists(fl):
                raise FileNotFoundError(f"{fl} not found: no f2 for {a} and {b}")
            m = read_rds(fl)
            vals = _vector(m)
            if len(vals) != 2 * nb:
                raise ValueError(f"{fl}: expected {nb} blocks")
            f2 = [math.nan if v is None else v for v in vals[:nb]]
            cnt = [0.0 if v is None else v * n for v, n in zip(vals[nb:], lengths)]  # noqa: B905
            means[i][j] = means[j][i] = f2
            counts[i][j] = counts[j][i] = cnt
    keep = list(range(nb))
    if remove_na:
        keep = [b for b in range(nb) if all(math.isfinite(means[i][j][b]) for i in range(len(pops))
                                            for j in range(i, len(pops)))]
        if not keep:
            raise ValueError("no blocks remain after discarding blocks with missing values")
    sel = [[[row[k] for k in keep] for row in rows] for rows in means]
    cnt = [[[row[k] for k in keep] for row in rows] for rows in counts]
    return _core.PyF2Blocks.from_block_means(pops, [lengths[k] for k in keep], sel, cnt)


def write_admixtools2_f2(f2_blocks: _core.PyF2Blocks, outdir: str | os.PathLike, *, overwrite: bool = False) -> None:
    """Write f2 blocks as an admixtools2 f2 directory (admixtools2 ``write_f2``),
    which admixtools2 reads with ``f2_from_precomp`` or ``read_f2``."""
    outdir = os.fspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    lengths_file = os.path.join(outdir, "block_lengths_f2.rds")
    if os.path.exists(lengths_file) and not overwrite:
        raise FileExistsError(f"{lengths_file} exists; pass overwrite=True to replace the directory's files")
    pops = list(f2_blocks.populations)
    write_rds(lengths_file, [int(x) for x in f2_blocks.block_lengths])
    for i, a in enumerate(pops):
        for b in pops[i:]:
            p1, p2 = _c_order(a, b)
            os.makedirs(os.path.join(outdir, p1), exist_ok=True)
            f2 = [float(x) for x in f2_blocks.f2_per_block(p1, p2)]
            lengths = list(f2_blocks.block_lengths)
            cnt = [c / n for c, n in zip(f2_blocks.counts_per_block(p1, p2), lengths)]  # noqa: B905
            write_rds(os.path.join(outdir, p1, f"{p2}_f2.rds"), matrix([f2, cnt], ["f2", "counts"]))


def extract_f2_subset(f2_blocks: _core.PyF2Blocks, populations: Sequence[str]) -> _core.PyF2Blocks:
    """These populations' f2 blocks (admixtools2 ``extract_f2_subset``)."""
    return f2_blocks.subset(list(populations))
