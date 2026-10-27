"""Search for the admixture graphs that best fit the data.

The searches are popmix's own, in Rust (``popmix._core.find_graphs``): every
graph with the given leaves and admixture nodes (``"exhaustive"``), or a
beam search that adds one leaf, then one admixture edge, at a time
(``"stepwise"``). Each candidate is fitted with qpgraph.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from popmix import _core

METHODS = ("stepwise", "exhaustive")


def find_graphs(f2_blocks: _core.F2Blocks, populations: Sequence[str], *, nadmix: int = 0,
                outgroup: str | None = None, method: str = "stepwise", **kwargs: Any):
    """Search for the admixture graphs that best fit the data.

    Parameters
    ----------
    f2_blocks
        f2 blocks with every population.
    populations
        The graph's leaves; ``populations[0]`` is the f3 base of every fit.
    nadmix
        Number of admixture nodes.
    outgroup
        Only graphs with this population directly below the root.
    method
        ``"stepwise"`` (default) or ``"exhaustive"``.
    **kwargs
        Options of the method (``beam_width`` for stepwise) and of qpgraph
        (``numstart``, ``diag``, ``diag_f3``, ``compat``).

    Returns
    -------
    tuple
        (graphs, n_evaluated): the graphs, best first, and the number of
        graphs fitted.
    """
    from popmix.errors import InputError
    from popmix.results import scored_graph

    if method not in METHODS:
        raise InputError(f"method must be one of {', '.join(METHODS)}, not {method!r}")
    graphs, n_evaluated = _core.find_graphs(f2_blocks, list(populations), nadmix, outgroup, method, **kwargs)
    return [scored_graph(g) for g in graphs], n_evaluated
