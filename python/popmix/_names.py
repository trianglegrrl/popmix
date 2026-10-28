"""Population arguments: a string is one population, as in admixtools2."""

from __future__ import annotations

from collections.abc import Iterable


def names(pops: str | Iterable[str]) -> list[str]:
    """``pops`` as a list of population names.

    Parameters
    ----------
    pops : str or iterable of str
        One population, or several.

    Returns
    -------
    list of str
        ``[pops]`` for a string (``list`` would split it into characters),
        else ``list(pops)``.
    """
    return [pops] if isinstance(pops, str) else list(pops)
