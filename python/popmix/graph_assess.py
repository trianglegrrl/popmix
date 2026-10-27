"""What a qpGraph fit can and cannot show (TASKS.md G9).

Every parameter added to a graph lowers its score, so a better score is not
evidence of a better graph (Maier et al. 2023, eLife, "on the limits of
fitting complex models of population history"). This module gives a fit's
out-of-sample score, fitting on alternate jackknife blocks and scoring on
the others, and warns when starts that reach the same best score disagree
on the admixture weights, which the f3 statistics then do not determine.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

from popmix.errors import PopmixWarning

if TYPE_CHECKING:
    from popmix import _core
    from popmix.results import QpGraphResult

# starts within this relative score of the best one count as equally good
SAME_SCORE = 1e-4
# a spread of the weights larger than this among them is not identifiable
WEIGHT_SPREAD = 0.05


def qpgraph_out_of_sample(f2_blocks: _core.F2Blocks, edges: Sequence[tuple[str, str]], **kwargs: Any) -> float:
    """Out-of-sample score of a graph: fitted on alternate jackknife blocks.

    The graph is fitted on the even-numbered blocks and scored on the
    odd-numbered ones (admixtools2's ``f2_blocks_test``), then the other way
    round; the result is the mean of the two scores. Unlike the score, it
    does not fall just because the graph has more parameters.

    Parameters
    ----------
    f2_blocks
        f2 blocks holding every leaf.
    edges
        The graph, as for ``qpgraph``.
    **kwargs
        Options of ``qpgraph``.

    Returns
    -------
    float
        The mean of the two out-of-sample scores.
    """
    from popmix.fits import qpgraph

    n = f2_blocks.n_blocks
    halves = [f2_blocks.select_blocks(list(range(start, n, 2))) for start in (0, 1)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PopmixWarning)
        scores = [qpgraph(train, edges, f2_blocks_test=test, **kwargs).score_test
                  for train, test in (halves, halves[::-1])]
    return sum(scores) / 2


def warn_unidentified(result: QpGraphResult) -> None:
    """Warn when equally good starts reach different admixture weights."""
    if len(result.opt) < 2 or not result.opt[0].weights:
        return
    best = min(o.score for o in result.opt)
    near = [o for o in result.opt if o.score <= best * (1 + SAME_SCORE) + 1e-12]
    for k in range(len(near[0].weights)):
        low, high = min(o.weights[k] for o in near), max(o.weights[k] for o in near)
        if high - low > WEIGHT_SPREAD:
            warnings.warn(f"qpgraph: {len(near)} starts reach the best score (within 0.01%) with weights from "
                          f"{low:.3f} to {high:.3f}: the admixture weights are not identifiable from these f3 "
                          "statistics", PopmixWarning, stacklevel=3)
            return
