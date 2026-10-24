"""Plots of popmix results. Needs matplotlib: ``pip install popmix[plot]``."""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pandas as pd


def _pyplot():
    try:
        import matplotlib

        if not os.environ.get("DISPLAY") and not os.environ.get("MPLBACKEND"):
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError("plotting needs matplotlib: pip install popmix[plot]") from e
    return plt


def heatmap(
    matrix: pd.DataFrame,
    *,
    title: str | None = None,
    center: float | None = None,
    cmap: str | None = None,
    annotate: bool | None = None,
    fmt: str = ".3g",
    path: str | os.PathLike[str] | None = None,
    ax: Any = None,
):
    """Heatmap of a population matrix (``outgroup_f3_matrix``, ``dstat_matrix``).

    Parameters
    ----------
    matrix
        Values to plot, with population labels as index and columns.
    title
        Plot title.
    center
        Value at the middle of a diverging color map (such as 0 for D or z);
        if None, a sequential map spans the data.
    cmap
        Matplotlib color map; default "RdBu_r" when centered, else "viridis".
    annotate
        Write each value in its cell; default when the matrix is at most 12
        populations across.
    fmt
        Format of the annotations.
    path
        If given, save the figure there (format from the extension).
    ax
        Axes to draw on; a new figure if None.

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _pyplot()
    values = matrix.to_numpy(dtype=float)
    if ax is None:
        n = max(matrix.shape)
        _, ax = plt.subplots(figsize=(2 + 0.5 * n, 1.5 + 0.5 * n))
    finite = values[np.isfinite(values)]
    if center is not None:
        span = np.max(np.abs(finite - center)) if finite.size else 1.0
        vmin, vmax = center - span, center + span
        cmap = cmap or "RdBu_r"
    else:
        vmin, vmax = (finite.min(), finite.max()) if finite.size else (0.0, 1.0)
        cmap = cmap or "viridis"
    image = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(matrix.shape[1]), labels=[str(c) for c in matrix.columns], rotation=90)
    ax.set_yticks(range(matrix.shape[0]), labels=[str(r) for r in matrix.index])
    ax.figure.colorbar(image, ax=ax)
    if annotate is None:
        annotate = max(matrix.shape) <= 12
    if annotate:
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                if np.isfinite(values[i, j]):
                    ax.text(j, i, format(values[i, j], fmt), ha="center", va="center", fontsize=7)
    if title:
        ax.set_title(title)
    ax.figure.tight_layout()
    if path is not None:
        ax.figure.savefig(path, dpi=150)
    return ax
