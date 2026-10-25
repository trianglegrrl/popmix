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


def pca(
    result,
    *,
    x: int = 1,
    y: int = 2,
    title: str | None = None,
    path: str | os.PathLike[str] | None = None,
    ax: Any = None,
):
    """Scatter plot of two principal components of a ``smartpca`` result.

    PCA samples are drawn as filled circles and projected samples as open
    triangles, colored by population.

    Parameters
    ----------
    result
        A ``popmix.PcaResult``.
    x, y
        Components on the axes (from 1).
    title
        Plot title.
    path
        If given, save the figure there.
    ax
        Axes to draw on; a new figure if None.

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _pyplot()
    df = result.coords
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 6))
    pops = list(dict.fromkeys(df["pop"]))
    colors = plt.get_cmap("tab20" if len(pops) > 10 else "tab10")
    for i, p in enumerate(pops):
        for role, marker, fill in (("pca", "o", True), ("project", "^", False)):
            d = df[(df["pop"] == p) & (df["role"] == role)]
            if d.empty:
                continue
            color = colors(i % colors.N)
            ax.scatter(d[f"PC{x}"], d[f"PC{y}"], marker=marker, s=18,
                       facecolors=color if fill else "none", edgecolors=color, label=p)
    ev = result.eigenvalues
    share = ev / ev.sum()
    ax.set_xlabel(f"PC{x} ({100 * share[x - 1]:.1f}%)")
    ax.set_ylabel(f"PC{y} ({100 * share[y - 1]:.1f}%)")
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))  # noqa: B905 (strict= needs Python 3.10)
    ax.legend(unique.values(), unique.keys(), fontsize=7, frameon=False, loc="best")
    if title:
        ax.set_title(title)
    ax.figure.tight_layout()
    if path is not None:
        ax.figure.savefig(path, dpi=150)
    return ax


def admixture(
    result,
    *,
    title: str | None = None,
    sort: bool = True,
    path: str | os.PathLike[str] | None = None,
    ax: Any = None,
):
    """Stacked bar plot of the ancestry proportions of an ``admixture`` result.

    Samples are grouped by population (in order of first appearance) and,
    with ``sort``, ordered within each by their largest component.

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _pyplot()
    df = result.Q
    k = result.k
    cols = [f"Q{j + 1}" for j in range(k)]
    pops = list(dict.fromkeys(df["pop"]))
    groups = []
    for p in pops:
        d = df[df["pop"] == p]
        if sort:
            main = d[cols].mean().idxmax()
            d = d.sort_values(main, ascending=False)
        groups.append(d)
    ordered = pd.concat(groups)
    if ax is None:
        _, ax = plt.subplots(figsize=(max(6, 0.08 * len(ordered)), 3))
    colors = plt.get_cmap("tab10" if k <= 10 else "tab20")
    bottom = np.zeros(len(ordered))
    x = np.arange(len(ordered))
    for j, c in enumerate(cols):
        ax.bar(x, ordered[c], bottom=bottom, width=1.0, color=colors(j % colors.N), linewidth=0)
        bottom += ordered[c].to_numpy()
    edges = np.cumsum([len(g) for g in groups])
    for e in edges[:-1]:
        ax.axvline(e - 0.5, color="black", lw=0.8)
    centers = edges - np.array([len(g) for g in groups]) / 2 - 0.5
    ax.set_xticks(centers, labels=pops, rotation=90)
    ax.set_xlim(-0.5, len(ordered) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_ylabel(f"K = {k}")
    if title:
        ax.set_title(title)
    ax.figure.tight_layout()
    if path is not None:
        ax.figure.savefig(path, dpi=150)
    return ax


def alder_curve(
    result,
    *,
    title: str | None = None,
    path: str | os.PathLike[str] | None = None,
    ax: Any = None,
):
    """The weighted LD curve of an ``alder`` result and its fitted exponential.

    Bins before the fit start are drawn lighter; the dashed line is the
    inter-chromosome affine term.

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _pyplot()
    c = result.curve
    df = c.curve[np.isfinite(c.curve["d_cm"])]
    fit = c.fit
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))
    before = df["d_cm"] < 100 * c.fit_start
    ax.scatter(df["d_cm"][before], df["wld"][before], s=4, color="0.7")
    ax.scatter(df["d_cm"][~before], df["wld"][~before], s=4, color="C0", label="weighted LD")
    d = df["d_cm"][~before].to_numpy()
    ax.plot(d, fit["amp_exp"] * np.exp(-fit["decay"] * d / 100) + fit["amp_aff"], color="C3",
            label=f"fit: {fit['decay']:.1f} ± {fit['decay_se']:.1f} generations")
    ax.axhline(fit["amp_aff"], color="C3", ls="--", lw=0.8)
    ax.set_xlabel("genetic distance (cM)")
    ax.set_ylabel("weighted LD")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(title or f"{result.test_pop} with {' and '.join(c.refs)}")
    ax.figure.tight_layout()
    if path is not None:
        ax.figure.savefig(path, dpi=150)
    return ax
