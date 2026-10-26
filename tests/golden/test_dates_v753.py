"""Golden master tests: popmix DATES vs DATES v753 run on the same data.

Reference outputs live in tests/golden/dates/v753/<config>/ and are produced by
generate_dates_reference.sh from the real DATES v753 binary. Each config holds
the DATES parameter file, so popmix is run with exactly the same settings.

Tolerances come from DATES itself, not from how close popmix happens to be:
- LD table (.out): pair counts exact; covariance, yreg and correlation to
  their printed precision (%15.9f, %12.6f, %12.6f).
- Fit estimates: printed precision (%9.3f) plus DATES's optimizer tolerance.
  DATES stops its GSL simplex within gslprecision = 1e-6 of r, which is up to
  1e-6 / step generations from the minimum; popmix minimises the same error
  to full precision. (chr20_22_q0: DATES stops at 3.164539 and prints 3.165;
  the true minimum, with lower error, is 3.164450.)
- Jackknife result: DATES computes it from per-chromosome estimates that are
  rounded and carry the optimizer error above, so the bound is the propagated
  error of those inputs (see jackknife_tol).

popmix corrects a few DATES quirks by default (DEVIATIONS.md, DATES-0001 to
DATES-0003). Every stage is checked with compat="dates-v753", which turns the
corrections off, at the tolerances above. The corrected default is then
checked three ways:
- default_ld_table: the corrections change only the bins they are meant to,
  and every other bin still matches DATES.
- default_result: where no correction changes the LD table, the result
  matches what DATES reports.
- default_vs_exact_mode: where a correction does change it, the corrected
  mesh agrees with DATES's exact mode, an independent oracle, to within the
  mesh's own approximation.
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pytest

import popmix as pm

REF_ROOT = Path(__file__).parent / "dates" / "v753"

# Half a unit in the last printed digit, plus a little room for binary representation.
COV_TOL = 5e-10 + 1e-15
CORR_TOL = 5e-7 + 1e-12
GEN_TOL = 5e-4 + 1e-9
DATES_GSL_PRECISION = 1e-6  # dates_expfit gslprecision, in r


def fit_tol(ref: "DatesReference") -> float:
    return GEN_TOL + DATES_GSL_PRECISION / ref.binsize


@dataclass
class DatesReference:
    name: str
    prefix: str
    source1: str
    source2: str
    test_pop: str
    binsize: float
    maxdis: float
    qbin: int
    lovalfit: float
    ld_rows: list[tuple[float, float, float, float, int]]  # (dist_cm, cov, yreg, corr, n_pairs)
    full_estimate: float
    jackknife_by_chrom: dict[int, tuple[float, float]]  # chrom -> (weight, estimate)
    jackknife_estimate: float
    jackknife_se: float


def _parse_par(path: Path) -> dict[str, str]:
    params = {}
    for line in path.read_text().splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            params[key.strip()] = value.strip()
    return params


def load_reference(name: str) -> DatesReference:
    ref_dir = REF_ROOT / name
    par = _parse_par(ref_dir / "par")
    source1, source2, test_pop, _ = (ref_dir / "admixlist").read_text().split()

    ld_rows = []
    for line in (ref_dir / "dates.out").read_text().splitlines():
        if line.lstrip().startswith("#") or not line.strip():
            continue
        dist, cov, yreg, corr, n_pairs = line.split()
        ld_rows.append((float(dist), float(cov), float(yreg), float(corr), int(n_pairs)))

    full_estimate = None
    for line in (ref_dir / "dates.flog").read_text().splitlines():
        if line.startswith("mean (generations):"):
            full_estimate = float(line.split(":")[1])
    assert full_estimate is not None, f"no full-data estimate in {name}/dates.flog"

    jackknife_by_chrom = {}
    for line in (ref_dir / "dates.jin").read_text().splitlines():
        chrom, weight, estimate = line.split()
        if float(weight) > 0:
            jackknife_by_chrom[int(chrom)] = (float(weight), float(estimate))

    jk_est, jk_se = (float(x) for x in (ref_dir / "dates.jout").read_text().split())

    # inputs kept in tests/golden/dates are named relative to the config directory
    prefix = str((ref_dir / par["genotypename"]).resolve())[: -len(".geno")]
    return DatesReference(
        name=name,
        prefix=prefix,
        source1=source1,
        source2=source2,
        test_pop=test_pop,
        binsize=float(par["binsize"]),
        maxdis=float(par["maxdis"]),
        qbin=int(par["qbin"]),
        lovalfit=float(par["lovalfit"]),
        ld_rows=ld_rows,
        full_estimate=full_estimate,
        jackknife_by_chrom=jackknife_by_chrom,
        jackknife_estimate=jk_est,
        jackknife_se=jk_se,
    )


def jackknife_tol(ref: DatesReference) -> float:
    """Bound on DATES's jackknife error from the error in its inputs.

    jackest = sum_k (mean - m_k) + sum_k w_k m_k / N. The mean and each m_k
    carry up to fit_tol of error, which adds up to at most (g + 2) * fit_tol
    over the g blocks, plus the final %9.3f rounding.
    """
    g = len(ref.jackknife_by_chrom)
    return (g + 2) * fit_tol(ref) + GEN_TOL


# yamnaya_chr21_22_q10: a source mixing pseudo-haploid samples and a diploid
# genome (DATES-0004), from AADR v66.p1 (generate_yamnaya_subset.sh)
QUICK_CONFIGS = ["chr20_22_q0", "chr20_22_q10", "yamnaya_chr21_22_q10"]
SLOW_CONFIGS = ["simulation_q0", "simulation_q10", "aadr_french_q10"]

COMPAT = "dates-v753"

# Mesh-mode configs whose exact-mode (qbin: 0) DATES run on the same data is
# an oracle for what the mesh approximates.
# Configs with a source of mixed ploidy, where DATES-0004 applies
MIXED_PLOIDY = {"yamnaya_chr21_22_q10"}

EXACT_MODE_ORACLE = {
    "chr20_22_q10": "chr20_22_q0",
    "simulation_q10": "simulation_q0",
}


_genotypes = {}


def _run_popmix(ref: DatesReference, compat: Optional[str]):
    if not os.path.exists(f"{ref.prefix}.geno"):
        pytest.skip(f"genotype data not available: {ref.prefix}")
    if ref.prefix not in _genotypes:
        _genotypes.clear()
        _genotypes[ref.prefix] = pm.read_genotypes(ref.prefix)
    return pm.dates(
        _genotypes[ref.prefix],
        test_pop=ref.test_pop,
        source1=ref.source1,
        source2=ref.source2,
        bin_width_cm=ref.binsize * 100.0,
        max_distance_cm=ref.maxdis * 100.0,
        min_distance_cm=ref.lovalfit,
        qbin=ref.qbin,
        compat=compat,
    )


_results = {}


def _result_for(name: str, compat: Optional[str] = COMPAT):
    """popmix result for a config.

    compat="dates-v753" reproduces DATES v753 exactly, quirks included, so
    every stage can be compared at DATES's own precision. compat=None is the
    corrected default, which is what users get.
    """
    key = (name, compat)
    if key not in _results:
        ref = load_reference(name)
        _results[key] = (ref, _run_popmix(ref, compat))
    return _results[key]


def _ld_mismatches(ref: DatesReference, ld_bins, skip=()) -> list[str]:
    bins = {round(b.distance_cm, 3): b for b in ld_bins}
    mismatches = []
    for k, (dist, cov, yreg, corr, n_pairs) in enumerate(ref.ld_rows):
        if k in skip:
            continue
        b = bins.get(round(dist, 3))
        if b is None:
            mismatches.append(f"{dist:.3f} cM: missing bin")
            continue
        if b.n_pairs != n_pairs:
            mismatches.append(f"{dist:.3f} cM: n_pairs {b.n_pairs} != {n_pairs}")
        if abs(b.weighted_ld - cov) > COV_TOL:
            mismatches.append(f"{dist:.3f} cM: cov {b.weighted_ld:.9f} != {cov:.9f}")
        if abs(b.regression_coef - yreg) > CORR_TOL:
            mismatches.append(f"{dist:.3f} cM: yreg {b.regression_coef:.6f} != {yreg:.6f}")
        if abs(b.correlation - corr) > CORR_TOL:
            mismatches.append(f"{dist:.3f} cM: corr {b.correlation:.6f} != {corr:.6f}")
    return mismatches


def _check_ld_table(name: str):
    ref, result = _result_for(name)
    mismatches = _ld_mismatches(ref, result.ld_bins)
    assert not mismatches, (
        f"{len(mismatches)} LD mismatches vs DATES v753 ({name}), first 10:\n"
        + "\n".join(mismatches[:10])
    )


def _check_full_fit(name: str):
    ref, result = _result_for(name)
    assert abs(result.fit.decay_rate - ref.full_estimate) <= fit_tol(ref), (
        f"full-data fit {result.fit.decay_rate:.6f} != DATES {ref.full_estimate:.3f}"
    )


def _check_jackknife_by_chrom(name: str):
    ref, result = _result_for(name)
    ours = {chrom: (weight, est) for chrom, weight, est in result.jackknife_by_chromosome}

    assert set(ours) == set(ref.jackknife_by_chrom), (
        f"chromosomes {sorted(ours)} != DATES {sorted(ref.jackknife_by_chrom)}"
    )
    mismatches = []
    for chrom, (weight, est) in ref.jackknife_by_chrom.items():
        our_weight, our_est = ours[chrom]
        if our_weight != weight:
            mismatches.append(f"chr{chrom}: weight {our_weight} != {weight}")
        if abs(our_est - est) > fit_tol(ref):
            mismatches.append(f"chr{chrom}: estimate {our_est:.6f} != {est:.3f}")
    assert not mismatches, "\n".join(mismatches)


def _check_jackknife_result(name: str, compat: Optional[str] = COMPAT):
    ref, result = _result_for(name, compat)
    tol = jackknife_tol(ref)
    assert abs(result.generations - ref.jackknife_estimate) <= tol, (
        f"jackknife estimate {result.generations:.6f} != DATES {ref.jackknife_estimate:.3f}"
        f" (tol {tol:.4f})"
    )
    assert abs(result.generations_se - ref.jackknife_se) <= tol, (
        f"jackknife SE {result.generations_se:.6f} != DATES {ref.jackknife_se:.3f}"
        f" (tol {tol:.4f})"
    )


def _corrected_bins(name: str) -> list[int]:
    """LD bins where the corrected default differs from compat="dates-v753"."""
    _, compat = _result_for(name, COMPAT)
    _, default = _result_for(name, None)
    assert len(compat.ld_bins) == len(default.ld_bins)
    return [
        k
        for k, (a, b) in enumerate(zip(compat.ld_bins, default.ld_bins))
        if (a.n_pairs, a.weighted_ld, a.correlation) != (b.n_pairs, b.weighted_ld, b.correlation)
    ]


def _check_default_ld_table(name: str):
    """The corrections change only the bins they are meant to.

    DATES-0001 and DATES-0002 are mesh-mode FFT quirks at lags past a
    chromosome's FFT window, so they can only change a tail of the table,
    from the first wrapped lag to maxdis, and nothing in exact mode. Every
    other bin of the default table must match DATES.
    """
    if name in MIXED_PLOIDY:
        pytest.skip("DATES-0004 changes the source frequencies, so every bin; see "
                    "every_sample_diploid_flag_counts_pseudohaploid_calls_as_two_alleles")
    ref, default = _result_for(name, None)
    corrected = _corrected_bins(name)
    n = len(default.ld_bins)
    if ref.qbin == 0:
        assert corrected == [], f"exact mode changed at bins {corrected[:10]}"
    else:
        assert corrected == list(range(n - len(corrected), n)), (
            f"corrected bins are not a tail of the table: {corrected[:10]}"
        )
    mismatches = _ld_mismatches(ref, default.ld_bins, skip=set(corrected))
    assert not mismatches, (
        f"{len(mismatches)} uncorrected bins differ from DATES v753 ({name}), first 10:\n"
        + "\n".join(mismatches[:10])
    )


def _check_default_result(name: str):
    """Where no correction changes the LD table, users get what DATES reports.

    DATES-0003 (fitting printed values) still applies, so this is the
    full-precision result against DATES's printed one.
    """
    if name in MIXED_PLOIDY:
        pytest.skip("DATES-0004 changes this config's source frequencies")
    if _corrected_bins(name):
        pytest.skip("DATES-0001 changes this config's LD table; see default_vs_exact_mode")
    _check_jackknife_result(name, compat=None)


def _check_default_vs_exact_mode(name: str):
    """The corrected mesh approximates DATES's exact mode, corrected bins included.

    The mesh bins a pair by its mesh-cell lag rather than its distance, so a
    pair within one cell (1/qbin of a bin) of a bin edge can land in the
    neighbouring bin: a bin's pair count can be off by up to 1/qbin of it.
    DATES's own mesh breaks this bound on chr20_22 (up to 113% off past
    91.7 cM, from wrapped pairs).

    The first bin is exempt: the mesh drops lag 0, which also holds each SNP
    paired with itself, so pairs within one cell are never counted (as in
    DATES). It lies below the fit range.
    """
    if name not in EXACT_MODE_ORACLE:
        pytest.skip(f"no exact-mode DATES run for {name}")
    ref, default = _result_for(name, None)
    exact = load_reference(EXACT_MODE_ORACLE[name])
    exact_pairs = {round(row[0], 3): row[4] for row in exact.ld_rows}
    bad = []
    for b in default.ld_bins[1:]:
        e = exact_pairs[round(b.distance_cm, 3)]
        if abs(b.n_pairs - e) > e / ref.qbin:
            bad.append(f"{b.distance_cm:.3f} cM: mesh {b.n_pairs} pairs vs exact {e}")
    assert not bad, f"{len(bad)} bins off by more than 1/qbin, first 10:\n" + "\n".join(bad[:10])


CHECKS = {
    "ld_table": _check_ld_table,
    "full_fit": _check_full_fit,
    "jackknife_by_chrom": _check_jackknife_by_chrom,
    "jackknife_result": _check_jackknife_result,
    "default_ld_table": _check_default_ld_table,
    "default_result": _check_default_result,
    "default_vs_exact_mode": _check_default_vs_exact_mode,
}


@pytest.mark.parametrize("check", list(CHECKS))
@pytest.mark.parametrize("name", QUICK_CONFIGS)
def test_dates_v753_quick(name: str, check: str):
    CHECKS[check](name)


@pytest.mark.slow
@pytest.mark.parametrize("check", list(CHECKS))
@pytest.mark.parametrize("name", SLOW_CONFIGS)
def test_dates_v753_live(name: str, check: str):
    if not (REF_ROOT / name / "dates.jout").exists():
        pytest.skip(f"DATES reference {name} not generated")
    CHECKS[check](name)
