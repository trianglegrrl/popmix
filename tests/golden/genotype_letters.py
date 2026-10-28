"""Genotypes as allele letters, for comparing files whose counted allele
may differ: a PED .map without alleles leaves each tool to choose it."""

from __future__ import annotations

import numpy as np


def letters(g) -> np.ndarray:
    """Each genotype as its two alleles, in sorted order ("" when missing)."""
    codes = g.to_numpy()
    first = np.array([s.allele_ref for s in g.snps])[:, None]
    second = np.array([s.allele_alt for s in g.snps])[:, None]
    swap = first > second
    lo, hi = np.where(swap, second, first), np.where(swap, first, second)
    het = np.char.add(lo, hi)
    return np.select([codes == 2, codes == 1, codes == 0],
                     [np.char.add(first, first), np.broadcast_to(het, codes.shape), np.char.add(second, second)], "")
