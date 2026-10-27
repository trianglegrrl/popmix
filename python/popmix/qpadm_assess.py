"""What a qpAdm result shows, from p, feasibility, the SEs and popdrop (TASKS.md G1).

p alone does not make a model: a model is only "not rejected", and only when
p exceeds alpha and every weight is in [0, 1] (Harney et al. 2021,
"Assessing the performance of qpAdm", Genetics 217). The notes flag what
weakens it: a weight within 2 SE of 0 (a source the data may not need) and
a model with fewer sources that is also not rejected (prefer the simplest).
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING

from popmix.errors import PopmixWarning

if TYPE_CHECKING:
    from popmix.results import PopDropRow, QpAdmResult


@dataclass(frozen=True)
class QpAdmAssessment:
    """What a qpAdm model's fit shows.

    ``status`` is one of:

    - ``"not tested"``: 0 degrees of freedom (as many right populations as
      left), so the data cannot reject the model;
    - ``"rejected"``: p at or below ``alpha``;
    - ``"infeasible"``: p above ``alpha``, but a weight outside [0, 1], which
      rejects the model too;
    - ``"not rejected"``: p above ``alpha`` and every weight in [0, 1]. Other
      models may fit as well.

    ``notes`` say what weakens a model that is not rejected (or infeasible).
    """

    status: str
    alpha: float
    notes: tuple[str, ...]


def assess(result: QpAdmResult, alpha: float = 0.05) -> QpAdmAssessment:
    """Assess a qpAdm result; see ``QpAdmResult.assess``."""
    if math.isnan(result.p):
        return QpAdmAssessment("not tested", alpha, ())
    if result.p <= alpha:
        return QpAdmAssessment("rejected", alpha, ())
    notes = []
    for source, w, se in zip(result.sources, result.weights, result.weight_se):
        if not 0.0 <= w <= 1.0:
            notes.append(f"{source}: weight {w:.3f} is outside [0, 1]")
        elif se > 0 and abs(w) < 2 * se:
            notes.append(f"{source}: weight {w:.3f} is within 2 SE ({se:.3f}) of 0")
    simpler = [row for row in result.popdrop if row.wt > 0 and row.feasible and row.p > alpha]
    if simpler:
        best = max(simpler, key=lambda row: (row.wt, row.p))
        model = f"{_sources(result.sources, best)} (p = {best.p:.3g})"
        if len(simpler) == 1:
            notes.append(f"a model with fewer sources is also not rejected: {model}")
        else:
            notes.append(f"{len(simpler)} models with fewer sources are also not rejected; the simplest: {model}")
    return QpAdmAssessment("not rejected" if result.feasible else "infeasible", alpha, tuple(notes))


def _sources(sources: tuple[str, ...], row: PopDropRow) -> str:
    kept = [s for s, flag in zip(sources, row.pattern) if flag == "0"]
    return f"{kept[0]} alone" if len(kept) == 1 else " + ".join(kept)


def warn_multiple_testing(n_models: int) -> None:
    """Say that a set of models' p-values are not corrected for testing them all."""
    if n_models > 1:
        warnings.warn(f"{n_models} qpAdm models tested: their p-values are not corrected for multiple testing, "
                      "so report how many models were tried with any model chosen (Harney et al. 2021)",
                      PopmixWarning, stacklevel=3)
