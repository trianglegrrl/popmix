"""What popmix's analyses return: frozen dataclasses.

They compare by value, pickle, and cannot be changed after the fact (their
lists are tuples, their arrays read-only). Statistics carry admixtools2's
names: ``est``, ``se``, ``z``, ``p``, ``chisq``, ``dof``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from popmix.qpadm_assess import QpAdmAssessment


def _array(rows: Any) -> np.ndarray:
    a = np.array(rows, dtype=float)
    a.flags.writeable = False
    return a


# --- qpWave and qpAdm ---------------------------------------------------------

@dataclass(frozen=True)
class RankTest:
    """One rank of a qpWave or qpAdm rank test."""

    rank: int
    chisq: float
    dof: int
    p: float


@dataclass(frozen=True)
class QpWaveResult:
    """qpWave: rank tests of the f4 matrix of ``left`` against ``right``.

    ``rank_tests`` go from rank 0 up; ``estimated_rank`` is the lowest rank
    not rejected at ``alpha`` (one above the highest rank tested if every
    test is rejected). Rank r means at least r + 1 streams of ancestry as
    seen from ``right``: a lower bound, since the rank is at most
    min(len(left), len(right)) - 1, right populations alike in their
    relation to ``left`` hide streams, and a test without power reads as
    "not rejected". Not rejecting rank 0 does not show a clade.
    """

    left: tuple[str, ...]
    right: tuple[str, ...]
    rank_tests: tuple[RankTest, ...]
    estimated_rank: int
    alpha: float
    f4_matrix: np.ndarray
    f4_covariance: np.ndarray

    def full_rank(self) -> RankTest | None:
        """The test of the highest rank.

        Returns
        -------
        RankTest or None
            The last of ``rank_tests``, or None if there are none.
        """
        return self.rank_tests[-1] if self.rank_tests else None


@dataclass(frozen=True)
class PopDropRow:
    """One model of qpAdm's ``popdrop`` table: the sources left out by
    ``pattern`` (admixtools2's columns)."""

    pattern: str
    wt: int
    dof: int
    chisq: float
    p: float
    f4rank: int
    weights: tuple[float, ...]
    feasible: bool
    best: bool | None
    dofdiff: int | None
    chisqdiff: float | None
    p_nested: float | None


@dataclass(frozen=True)
class QpAdmResult:
    """qpAdm: ``target`` as a mixture of ``sources``, tested against the
    right populations."""

    target: str
    sources: tuple[str, ...]
    weights: tuple[float, ...]
    weight_se: tuple[float, ...]
    chisq: float
    dof: int
    p: float
    constrained: bool
    popdrop: tuple[PopDropRow, ...]
    #: Rank tests of the f4 matrix of the target and sources against the
    #: right populations, every rank (admixtools2 ``rankdrop``).
    rankdrop: tuple[RankTest, ...] = ()
    #: Jackknife blocks and SNPs (the block lengths' sum) of the f2 blocks.
    n_blocks: int = 0
    n_snps: int = 0
    #: The right populations.
    right: tuple[str, ...] = ()

    @property
    def weight_z(self) -> tuple[float, ...]:
        """Each weight over its SE (NaN where the SE is 0)."""
        return tuple(w / se if se > 0 else float("nan") for w, se in zip(self.weights, self.weight_se))

    @property
    def feasible(self) -> bool:
        """Whether every weight is in [0, 1]."""
        return all(0.0 <= w <= 1.0 for w in self.weights)

    def weights_dict(self) -> dict[str, float]:
        """The weights by source, in source order.

        Returns
        -------
        dict of str to float
            Source population to its weight.
        """
        return dict(zip(self.sources, self.weights))

    def assess(self, alpha: float = 0.05) -> QpAdmAssessment:
        """What the fit shows, from p, feasibility, the SEs and popdrop.

        p alone does not make a model good: see ``QpAdmAssessment``.

        Parameters
        ----------
        alpha
            Significance level: p at or below it rejects the model.

        Returns
        -------
        QpAdmAssessment
            The status (not tested, rejected, infeasible or not rejected)
            and notes: weights within 2 SE of 0 or outside [0, 1], and
            models with fewer sources that are also not rejected.
        """
        from popmix.qpadm_assess import assess

        return assess(self, alpha)


@dataclass(frozen=True)
class F4RatioResult:
    """qpF4ratio: the share ``alpha`` of ``source_a`` in ``test``,
    f4(A, O; X, C) / f4(A, O; B, C) for test X, sources A and B, outgroup O
    and reference C, with its jackknife SE."""

    test: str
    source_a: str
    source_b: str
    outgroup: str
    reference: str
    alpha: float
    se: float
    z: float
    n_blocks: int
    f4_numerator: float
    f4_denominator: float


# --- qpGraph -------------------------------------------------------------------

@dataclass(frozen=True)
class QpGraphEdge:
    """A fitted edge: a drift length, or an admixture weight (``kind``)."""

    from_node: str
    to_node: str
    kind: str
    weight: float
    low: float | None
    high: float | None


@dataclass(frozen=True)
class QpGraphF3:
    """An f3 statistic, as estimated and as the graph fits it."""

    pop1: str
    pop2: str
    pop3: str
    est: float
    se: float
    fit: float
    diff: float
    z: float
    p: float


@dataclass(frozen=True)
class QpGraphStart:
    """One optimization start of the admixture weights (admixtools2 ``opt``)."""

    start: tuple[float, ...]
    weights: tuple[float, ...]
    score: float
    evaluations: int
    convergence: int

    @property
    def converged(self) -> bool:
        return self.convergence == 0


@dataclass(frozen=True)
class QpGraphResult:
    """qpGraph: an admixture graph fitted to the f2 blocks."""

    pops: tuple[str, ...]
    score: float
    edges: tuple[QpGraphEdge, ...]
    f3: tuple[QpGraphF3, ...]
    opt: tuple[QpGraphStart, ...]
    ppinv: np.ndarray
    score_test: float | None

    @property
    def n_f3(self) -> int:
        """The f3 statistics fitted: f3(base; i, j) for i <= j."""
        return len(self.f3)

    @property
    def n_params(self) -> int:
        """Free parameters: the drift lengths and the admixture weights."""
        drift = sum(e.kind != "admix" for e in self.edges)
        return drift + sum(e.kind == "admix" for e in self.edges) // 2

    @property
    def dof(self) -> int:
        """f3 statistics less parameters; at or below 0 the graph can fit anything."""
        return self.n_f3 - self.n_params

    @property
    def penalized_score(self) -> float:
        """The score plus twice the parameters (AIC-style), for comparing graphs
        with different numbers of admixtures: the score falls with every one
        added."""
        return self.score + 2 * self.n_params


@dataclass(frozen=True)
class ScoredGraph:
    """A graph found by ``find_graphs``, with its qpgraph score. The path from
    the root to the first population comes first in ``edges``: it is the f3
    base."""

    edges: tuple[tuple[str, str], ...]
    score: float


# --- DATES ---------------------------------------------------------------------

@dataclass(frozen=True)
class LdBin:
    """Weighted LD in one distance bin."""

    distance_cm: float
    n_pairs: int
    weighted_ld: float
    se: float
    z: float
    correlation: float
    regression_coef: float


@dataclass(frozen=True)
class DecayFit:
    """The exponential fitted to the LD decay."""

    amplitude: float
    decay_rate: float
    decay_rate_se: float
    intercept: float
    r_squared: float
    rss: float


@dataclass(frozen=True)
class PulseCheck:
    """Whether two admixture pulses fit the LD decay better than one.

    A₁·e^(−n₁·d) + A₂·e^(−n₂·d) + K fitted to the same bins as the single
    pulse, with leave-one-chromosome-out jackknife SEs. ``multiple_pulses``
    is True only when both amplitudes and n₂ − n₁ are more than 2 SEs from
    zero: then the single date is an average over admixture at different
    times. False means the curve gives no evidence against one pulse, which
    is not evidence for one: admixture over a short span fits one
    exponential too.
    """

    generations: tuple[float, float]
    generations_se: tuple[float, float]
    amplitudes: tuple[float, float]
    amplitudes_se: tuple[float, float]
    difference_se: float
    rss: float
    multiple_pulses: bool


@dataclass(frozen=True)
class DatesResult:
    """DATES: the date of admixture of ``admixed`` from its two ``sources``.

    ``generations`` and ``generations_se`` are the weighted jackknife
    estimate and SE that DATES reports; ``fit.decay_rate`` is the full-data
    fit; ``jackknife_by_chromosome`` holds (chromosome, estimate, weight).
    ``years_bp`` is ``sample_age`` + generations x ``generation_time``:
    years before present when the samples' age is given, else before they
    lived; ``years_bp_se`` combines the SEs of all three.
    ``source_alleles``: each source's mean observed alleles per SNP used
    (a warning below 10).
    """

    admixed: str
    sources: tuple[str, ...]
    generations: float
    generations_se: float
    years_bp: float
    years_bp_se: float
    generation_time: float
    generation_time_se: float
    sample_age: float
    sample_age_se: float
    n_snps: int
    n_chromosomes: int
    source_alleles: tuple[float, float]
    jackknife_by_chromosome: tuple[tuple[int, float, float], ...]
    ld_bins: tuple[LdBin, ...]
    fit: DecayFit
    #: The two-pulse check (``PulseCheck``), None without the jackknife.
    pulse_check: PulseCheck | None = None

    def to_ld_dataframe(self):
        """The LD bins as a pandas DataFrame, one row per bin.

        Returns
        -------
        pandas.DataFrame
            ``distance_cm``, ``n_pairs``, ``weighted_ld``, ``se``, ``z``,
            ``correlation`` and ``regression_coef`` (the ``LdBin`` fields),
            one row per bin in ``ld_bins`` order.
        """
        import dataclasses

        import pandas as pd

        columns = [f.name for f in dataclasses.fields(LdBin)]
        return pd.DataFrame([dataclasses.astuple(b) for b in self.ld_bins], columns=columns)


# --- from the Rust core's objects ---------------------------------------------

def _fields(cls, obj, **override) -> Any:
    import dataclasses

    values = {f.name: override[f.name] if f.name in override else getattr(obj, f.name)
              for f in dataclasses.fields(cls)}
    for name, value in values.items():
        if isinstance(value, list):
            values[name] = tuple(value)
    return cls(**values)


def qpwave_result(r) -> QpWaveResult:
    return _fields(QpWaveResult, r, rank_tests=tuple(_fields(RankTest, t) for t in r.rank_tests),
                   f4_matrix=_array(r.f4_matrix), f4_covariance=_array(r.f4_covariance))


def qpadm_result(r, f2_blocks=None, right=(), counts=(0, 0)) -> QpAdmResult:
    counts = {"n_blocks": counts[0], "n_snps": counts[1]} if f2_blocks is None else \
        {"n_blocks": f2_blocks.n_blocks, "n_snps": sum(f2_blocks.block_lengths)}
    return _fields(QpAdmResult, r, popdrop=tuple(_fields(PopDropRow, row) for row in r.popdrop),
                   rankdrop=tuple(_fields(RankTest, t) for t in r.rankdrop), right=tuple(right), **counts)


def f4ratio_result(r) -> F4RatioResult:
    return _fields(F4RatioResult, r)


def qpgraph_result(r) -> QpGraphResult:
    return _fields(QpGraphResult, r, edges=tuple(_fields(QpGraphEdge, e) for e in r.edges),
                   f3=tuple(_fields(QpGraphF3, f) for f in r.f3),
                   opt=tuple(_fields(QpGraphStart, s) for s in r.opt), ppinv=_array(r.ppinv))


def scored_graph(g) -> ScoredGraph:
    return ScoredGraph(edges=tuple(tuple(e) for e in g.edges), score=g.score)


def dates_result(r) -> DatesResult:
    check = r.pulse_check
    return _fields(DatesResult, r, ld_bins=tuple(_fields(LdBin, b) for b in r.ld_bins),
                   fit=_fields(DecayFit, r.fit),
                   pulse_check=None if check is None else _fields(PulseCheck, check),
                   jackknife_by_chromosome=tuple(tuple(j) for j in r.jackknife_by_chromosome))
