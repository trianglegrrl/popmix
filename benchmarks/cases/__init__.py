"""Benchmark cases: popmix and a reference tool on identical inputs.

Every case runs the popmix command line (one process, as a user would) and
the reference tool's own program, each writing its result to its own
directory; `agree` then compares the two results at the CLAUDE.md
tolerances. Inputs are made once, by the reference tools, in the cache.
"""

from __future__ import annotations

from ._common import AADR, CACHE, CONVERTF, PLINK2, POPS, SMARTPCA, _f2_mismatches, _r, subset
from .classic import (
    CLASSIC_QP3POP_CASE,
    CLASSIC_QPADM_CASE,
    CLASSIC_QPDSTAT_CASE,
    CLASSIC_QPGRAPH_CASE,
)
from .dating import DATES_CASE, ROLLOFF_CASE
from .fstats import (
    _R_F2,
    CONVERTF_AADR_CASE,
    CONVERTF_CASE,
    F2_CASE,
    QPADM_CASE,
    QPDSTAT_CASE,
    RESAMPLE_CASE,
    ROTATE_CASE,
    _agree_f2,
)
from .graphs import _R_FIND_GRAPHS, _R_QPGRAPH, FIND_GRAPHS_CASE, QPGRAPH_CASE
from .samples import (
    ADMIXTURE_CASE,
    KING_CASE,
    PCA_CASE,
    READ_CASE,
    SPREAD_SEEDS,
    _agree_king,
    _agree_pca,
    _seed_runs,
)

# the names used from outside the package (run.py, scaling.py, the tests)
__all__ = [
    "AADR",
    "CACHE",
    "CASES",
    "CONVERTF",
    "F2_CASE",
    "PLINK2",
    "POPS",
    "SMARTPCA",
    "SPREAD_SEEDS",
    "_R_F2",
    "_R_FIND_GRAPHS",
    "_R_QPGRAPH",
    "_agree_f2",
    "_agree_king",
    "_agree_pca",
    "_f2_mismatches",
    "_r",
    "_seed_runs",
    "subset",
]

CASES = {c.name: c for c in (CONVERTF_CASE, CONVERTF_AADR_CASE, F2_CASE, QPADM_CASE, PCA_CASE, KING_CASE, QPGRAPH_CASE,
                             QPDSTAT_CASE, ADMIXTURE_CASE, ROLLOFF_CASE, DATES_CASE, READ_CASE,
                             CLASSIC_QPDSTAT_CASE, CLASSIC_QP3POP_CASE, CLASSIC_QPADM_CASE, CLASSIC_QPGRAPH_CASE,
                             ROTATE_CASE, FIND_GRAPHS_CASE, RESAMPLE_CASE)}
