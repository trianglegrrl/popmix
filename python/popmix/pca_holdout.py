"""Shrinkage of projected samples, measured on held-out reference samples.

A PCA sample's coordinates include its own noise, which helped fit the axes;
a projected sample's do not, so projected samples sit closer to the origin
(EIGENSOFT's documentation gives the test: two halves of one population, one
in the PCA and one projected, should coincide). With ``holdout`` the halves
are made, and each component's shrinkage is the least-squares slope of the
held-out halves' population means on the PCA halves', both centered: a
measurement, not a model.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from popmix.errors import PopmixWarning

#: Shrinkage below this warns.
WARN_BELOW = 0.9


def measured_shrinkage(coords: pd.DataFrame) -> np.ndarray:
    """Each component's shrinkage from ``coords`` with "holdout" rows."""
    pcs = [c for c in coords.columns if c.startswith("PC")]
    means = coords[coords.role.isin(["pca", "holdout"])].groupby(["pop", "role"])[pcs].mean().unstack("role")
    means = means.dropna()
    a = means.xs("pca", axis=1, level="role").to_numpy()
    b = means.xs("holdout", axis=1, level="role").to_numpy()
    a, b = a - a.mean(0), b - b.mean(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (a * b).sum(0) / (a * a).sum(0)


def warn_shrinkage(shrinkage: np.ndarray) -> None:
    """Warn where the projected samples are pulled in by more than 10%."""
    low = [(j, s) for j, s in enumerate(shrinkage, 1) if s < WARN_BELOW]
    if low:
        warnings.warn("projected samples sit closer to the origin than the PCA samples: held-out reference samples "
                      "reach " + ", ".join(f"{s:.0%} (PC{j})" for j, s in low) + " of their populations' spread; "
                      "compare projected samples with the held-out ones (role \"holdout\"), not with the PCA samples",
                      PopmixWarning, stacklevel=3)
