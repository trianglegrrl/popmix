"""Models fitted to f2 blocks: qpWave, qpAdm, qpF4ratio and qpGraph."""

from __future__ import annotations

from collections.abc import Sequence

from popmix import _core
from popmix.results import (
    F4RatioResult,
    QpAdmResult,
    QpGraphResult,
    QpWaveResult,
    f4ratio_result,
    qpadm_result,
    qpgraph_result,
    qpwave_result,
)


def qpwave(f2_blocks: _core.F2Blocks, left: Sequence[str], right: Sequence[str], *,
           fudge: float = 1e-4, alpha: float = 0.05) -> QpWaveResult:
    """qpWave (admixtools2 ``qpwave``): how many waves of ancestry, as seen
    from ``right``, ``left`` descends from.

    Parameters
    ----------
    f2_blocks
        f2 blocks holding every population.
    left, right
        The left and right (outgroup) populations.
    fudge
        Added to the diagonal of the f4 covariance, times its trace.
    alpha
        Significance level of the rank tests: the estimated rank is the
        lowest rank whose ``p`` exceeds it.

    Returns
    -------
    QpWaveResult
        The rank tests, rank 0 first, and the estimated rank: a lower
        bound on the number of streams minus one (see ``QpWaveResult``).
    """
    return qpwave_result(_core.qpwave(f2_blocks, list(left), list(right), fudge=fudge, alpha=alpha))


def qpadm(f2_blocks: _core.F2Blocks, left: Sequence[str], right: Sequence[str], target: str, *,
          constrained: bool = False, fudge: float = 1e-4, compat: str | None = None) -> QpAdmResult:
    """qpAdm (admixtools2 ``qpadm``): ``target`` as a mixture of the ``left``
    populations, tested against ``right``.

    Parameters
    ----------
    f2_blocks
        f2 blocks holding every population.
    left
        The sources (admixtools2's ``left``, the target apart).
    right
        The right (outgroup) populations.
    target
        The population modelled.
    constrained
        Keep the weights in [0, 1].
    fudge
        Added to the diagonal of the f4 covariance, times its trace.
    compat
        ``"admixtools2-2.0.10"`` reproduces admixtools2's known quirks
        instead of correcting them (DEVIATIONS.md); None corrects them.

    Returns
    -------
    QpAdmResult
        Weights, their SEs, the model's chisq, dof and p, and popdrop.
    """
    return qpadm_result(_core.qpadm(f2_blocks, list(left), list(right), target,
                                    constrained=constrained, fudge=fudge, compat=compat))


def qpf4ratio(f2_blocks: _core.F2Blocks, test: str, source_a: str, source_b: str, outgroup: str,
              reference: str) -> F4RatioResult:
    """qpF4ratio: the share alpha of ``source_a`` in ``test``,
    f4(A, O; X, C) / f4(A, O; B, C) for test X, sources A and B, outgroup O
    and reference C.

    Parameters
    ----------
    f2_blocks
        f2 blocks holding all five populations.
    test
        The admixed population X.
    source_a
        The source A whose share is estimated.
    source_b
        The other source B.
    outgroup
        The outgroup O.
    reference
        The reference C.

    Returns
    -------
    F4RatioResult
        alpha, its jackknife SE and z, and the two f4 statistics.
    """
    return f4ratio_result(_core.qpf4ratio(f2_blocks, test, source_a, source_b, outgroup, reference))


def qpgraph(f2_blocks: _core.F2Blocks, edges: Sequence[tuple[str, str]], *,
            admix_weights: Sequence[float] | None = None, numstart: int = 10, seed: int = 1,
            diag: float = 1e-4, diag_f3: float = 1e-5, constrained: bool = True,
            f2_blocks_test: _core.F2Blocks | None = None, compat: str | None = None) -> QpGraphResult:
    """Fit an admixture graph, as admixtools2 ``qpgraph`` does.

    Drift lengths minimise the covariance-weighted residual of f3(base; i, j)
    for every pair of leaves, with lengths >= 0; admixture weights are
    optimized over [0, 1] from ``numstart`` random starts. The populations are
    the graph's leaves in order of first appearance in ``edges``; the first is
    the f3 base.

    Parameters
    ----------
    f2_blocks
        f2 blocks holding every leaf.
    edges
        The graph as (from, to) pairs; a node with two incoming edges is an
        admixture node.
    admix_weights
        Fit at these weights instead of optimizing them, one per admixture
        node in order of first appearance (the weight of the incoming edge
        whose source appears first).
    numstart
        Random starts for the admixture weights.
    seed
        Seed for the starts, as R's ``set.seed``: the same seed gives
        admixtools2's starts.
    diag
        Regularization of the edge-length fit.
    diag_f3
        Regularization of the f3 covariance.
    constrained
        Lengths >= 0 and weights in [0, 1].
    f2_blocks_test
        Also score the fit against these blocks' f3 statistics
        (``score_test``), for an out-of-sample score.
    compat
        None polishes the best start's weights past a zero-bound kink, where
        R's L-BFGS-B (and so admixtools2) stalls; "admixtools2-2.0.10" stops
        there, as admixtools2 (ADMIXTOOLS2-0009). The starts in ``opt`` are
        unpolished either way.

    Returns
    -------
    QpGraphResult
        Fitted edges, score, f3 fits and the optimization starts.
    """
    from popmix.graph_assess import warn_unidentified

    weights = None if admix_weights is None else list(admix_weights)
    result = qpgraph_result(_core.qpgraph(
        f2_blocks, [tuple(e) for e in edges], admix_weights=weights, numstart=numstart, seed=seed,
        diag=diag, diag_f3=diag_f3, constrained=constrained, f2_blocks_test=f2_blocks_test, compat=compat))
    warn_unidentified(result)
    return result
