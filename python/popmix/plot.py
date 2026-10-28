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

    PCA samples are drawn as filled circles, held-out reference samples as
    open circles and projected samples as open triangles, colored by
    population.

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
        for role, marker, fill in (("pca", "o", True), ("holdout", "o", False), ("project", "^", False)):
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


def rolloff_curve(
    result,
    *,
    title: str | None = None,
    path: str | os.PathLike[str] | None = None,
    ax: Any = None,
):
    """The weighted LD correlation of a ``rolloff`` result and its fitted
    exponential C + A exp(-m d).

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _pyplot()
    b = result.bins[np.isfinite(result.bins["corr"])]
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(b["d_cm"], b["corr"], s=4, color="C0", label="weighted LD correlation")
    d = b["d_cm"].to_numpy()
    mean, se = result.decay_jackknife
    ax.plot(d, result.affine + result.amplitude * np.exp(-result.decay * d / 100), color="C3",
            label=f"fit: {mean:.1f} ± {se:.1f} generations")
    ax.axhline(result.affine, color="C3", ls="--", lw=0.8)
    ax.set_xlabel("genetic distance (cM)")
    ax.set_ylabel("correlation")
    ax.legend(frameon=False, fontsize=8)
    if title:
        ax.set_title(title)
    ax.figure.tight_layout()
    if path is not None:
        ax.figure.savefig(path, dpi=150)
    return ax


def _graph_layout(edges: list[tuple[str, str]]) -> dict[str, tuple[float, float]]:
    """Node positions: y from 1 at the root down to 0, each node one step
    below its lowest parent (the longest path from the root) and the leaves
    at 0; x with the leaves in order along the first-parent tree and each
    internal node above its children."""
    children: dict[str, list[str]] = {}
    parents: dict[str, list[str]] = {}
    for a, b in edges:
        children.setdefault(a, []).append(b)
        parents.setdefault(b, []).append(a)
    root = next(a for a, _ in edges if a not in parents)
    depth: dict[str, int] = {root: 0}
    pending = [root]
    while pending:  # longest path, by relaxing until no depth grows
        v = pending.pop()
        for c in children.get(v, []):
            if depth.get(c, -1) < depth[v] + 1:
                depth[c] = depth[v] + 1
                pending.append(c)
    bottom = max(depth.values()) or 1
    y = {v: 0.0 if not children.get(v) else 1 - d / bottom for v, d in depth.items()}
    x: dict[str, float] = {}
    next_leaf = [0.0]

    def place(v: str) -> float:
        kids = [c for c in children.get(v, []) if parents[c][0] == v]
        if not children.get(v):
            x[v] = next_leaf[0]
            next_leaf[0] += 1
        elif kids:
            x[v] = float(np.mean([place(c) for c in kids]))
        else:  # only admixture children, placed from their first parent
            x[v] = next_leaf[0]
            next_leaf[0] += 1
        return x[v]

    place(root)
    return {v: (x[v], y[v]) for v in y}


def graph(
    graph_or_fit,
    *,
    title: str | None = None,
    path: str | os.PathLike[str] | None = None,
    ax: Any = None,
):
    """An admixture graph: a ``qpgraph`` result (edges labelled as
    admixtools2's plots are: drift length x 1000, admixture weight in %), or
    the graph as (from, to) pairs.

    Admixture edges are dashed; leaves are labelled with their names.

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _pyplot()
    if hasattr(graph_or_fit, "edges") and not isinstance(graph_or_fit, (list, tuple)):
        rows = [(e.from_node, e.to_node, e.kind, e.weight) for e in graph_or_fit.edges]
    else:
        rows = [(a, b, None, None) for a, b, *_ in graph_or_fit]
    edges = [(a, b) for a, b, _, _ in rows]
    indeg: dict[str, int] = {}
    for _, b in edges:
        indeg[b] = indeg.get(b, 0) + 1
    pos = _graph_layout(edges)
    leaves = {b for _, b in edges} - {a for a, _ in edges}
    if ax is None:
        _, ax = plt.subplots(figsize=(max(5, 0.9 * len(leaves)), 5))
    for a, b, kind, w in rows:
        admix = (kind == "admix") if kind is not None else indeg[b] > 1
        (x0, y0), (x1, y1) = pos[a], pos[b]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops={"arrowstyle": "-|>", "color": "C3" if admix else "0.3",
                                "linestyle": "--" if admix else "-", "lw": 1})
        if w is not None and np.isfinite(w):
            label = f"{round(w * 100)}%" if admix else f"{round(w * 1000)}"
            ax.text((x0 + x1) / 2, (y0 + y1) / 2, label, fontsize=7, ha="center", va="center",
                    color="C3" if admix else "0.2", bbox={"boxstyle": "round,pad=0.1", "fc": "white", "lw": 0})
    for v, (xv, yv) in pos.items():
        if v in leaves:
            ax.text(xv, yv - 0.15, v, ha="center", va="top", fontsize=8, rotation=30)
        else:
            ax.plot(xv, yv, "o", color="0.3", ms=3)
    ys = [p[1] for p in pos.values()]
    xs = [p[0] for p in pos.values()]
    ax.set_xlim(min(xs) - 0.8, max(xs) + 0.8)
    ax.set_ylim(min(ys) - 1.2, max(ys) + 0.5)
    ax.axis("off")
    score = getattr(graph_or_fit, "score", None)
    ax.set_title(title or (f"score {score:.3f}" if score is not None else ""))
    ax.figure.tight_layout()
    if path is not None:
        ax.figure.savefig(path, dpi=150)
    return ax


def qpgraph_fit(
    fit,
    *,
    title: str | None = None,
    path: str | os.PathLike[str] | None = None,
    ax: Any = None,
):
    """How well a ``qpgraph`` result fits each f3-statistic: its residual
    (estimate - fit) / SE, with |z| = 3 marked.

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _pyplot()
    rows = [(f"{f.pop1}; {f.pop2}, {f.pop3}", f.z) for f in fit.f3 if np.isfinite(f.z)]
    if ax is None:
        _, ax = plt.subplots(figsize=(6, max(2.5, 0.22 * len(rows))))
    yv = np.arange(len(rows))
    z = np.array([r[1] for r in rows])
    ax.scatter(z, yv, s=12, color=np.where(np.abs(z) > 3, "C3", "C0"))
    for v in (-3, 3):
        ax.axvline(v, color="0.5", ls="--", lw=0.8)
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_yticks(yv, labels=[r[0] for r in rows], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("f3 residual z")
    ax.set_title(title or f"score {fit.score:.3f}, worst |z| {np.max(np.abs(z)):.2f}")
    ax.figure.tight_layout()
    if path is not None:
        ax.figure.savefig(path, dpi=150)
    return ax
