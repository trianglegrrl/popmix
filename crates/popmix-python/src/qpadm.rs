//! qpWave, qpAdm and f4-ratio for Python: single models, rotations,
//! pairwise cladality and model sets.

use crate::compute;
use crate::f2::PyF2Blocks;
use popmix_core::python::input_error;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Result of a single rank test from qpWave.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRankTest {
    #[pyo3(get)]
    rank: usize,
    #[pyo3(get)]
    chi_sq: f64,
    #[pyo3(get)]
    df: usize,
    #[pyo3(get)]
    p_value: f64,
}

/// Result of qpWave analysis.
#[pyclass]
pub(crate) struct PyQpWaveResult {
    #[pyo3(get)]
    rank_tests: Vec<PyRankTest>,
    #[pyo3(get)]
    estimated_rank: usize,
    #[pyo3(get)]
    left_pops: Vec<String>,
    #[pyo3(get)]
    right_pops: Vec<String>,
    /// f4(left[i+1], left[0]; right[j+1], right[0]), as rows i, columns j.
    #[pyo3(get)]
    f4_matrix: Vec<Vec<f64>>,
    /// Jackknife covariance of the f4 matrix flattened by rows (left first).
    #[pyo3(get)]
    f4_covariance: Vec<Vec<f64>>,
}

#[pymethods]
impl PyQpWaveResult {
    /// Check if populations form a clade (rank 0 p-value > threshold).
    #[pyo3(signature = (p_threshold=0.05))]
    fn is_clade(&self, p_threshold: f64) -> bool {
        self.rank_tests
            .first()
            .map(|r| r.p_value > p_threshold)
            .unwrap_or(false)
    }

    /// Get the full rank test.
    fn full_rank(&self) -> Option<PyRankTest> {
        self.rank_tests.last().cloned()
    }
}

/// Result of f4 ratio analysis for admixture proportion estimation.
#[pyclass]
pub(crate) struct PyF4RatioResult {
    #[pyo3(get)]
    test: String,
    #[pyo3(get)]
    source_a: String,
    #[pyo3(get)]
    source_b: String,
    #[pyo3(get)]
    outgroup: String,
    #[pyo3(get)]
    reference: String,
    #[pyo3(get)]
    alpha: f64,
    #[pyo3(get)]
    std_err: f64,
    #[pyo3(get)]
    z_score: f64,
    #[pyo3(get)]
    n_blocks: usize,
    #[pyo3(get)]
    f4_numerator: f64,
    #[pyo3(get)]
    f4_denominator: f64,
}

#[pymethods]
impl PyF4RatioResult {
    fn __repr__(&self) -> String {
        format!(
            "F4RatioResult(test='{}', source_a='{}', source_b='{}', alpha={:.4}, std_err={:.4})",
            self.test, self.source_a, self.source_b, self.alpha, self.std_err
        )
    }
}

/// Estimate admixture proportions via f4 ratios.
///
/// Computes alpha = f4(A,O;X,C) / f4(A,O;B,C) where:
/// - X is the test (admixed) population
/// - A is source population contributing alpha proportion
/// - B is source population contributing (1-alpha) proportion
/// - O is an outgroup
/// - C is a reference population
///
/// Args:
///     f2_blocks: Pre-computed F2 blocks
///     test: Test (admixed) population X
///     source_a: Source population A
///     source_b: Source population B
///     outgroup: Outgroup population O
///     reference: Reference population C
///
/// Returns:
///     PyF4RatioResult with estimated alpha and standard error
#[pyfunction]
#[pyo3(signature = (f2_blocks, test, source_a, source_b, outgroup, reference))]
pub(crate) fn qpf4ratio(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    test: String,
    source_a: String,
    source_b: String,
    outgroup: String,
    reference: String,
) -> PyResult<PyF4RatioResult> {
    compute(py, move || {
        let result = popmix_core::f4ratio(
            &f2_blocks.inner,
            &test,
            &source_a,
            &source_b,
            &outgroup,
            &reference,
            None,
        )?;

        Ok(PyF4RatioResult {
            test: result.test.to_string(),
            source_a: result.source_a.to_string(),
            source_b: result.source_b.to_string(),
            outgroup: result.outgroup.to_string(),
            reference: result.reference.to_string(),
            alpha: result.alpha,
            std_err: result.std_err,
            z_score: result.z_score,
            n_blocks: result.n_blocks,
            f4_numerator: result.f4_numerator,
            f4_denominator: result.f4_denominator,
        })
    })
}

/// Result of qpAdm analysis.
#[pyclass]
pub(crate) struct PyQpAdmResult {
    #[pyo3(get)]
    target: String,
    #[pyo3(get)]
    sources: Vec<String>,
    #[pyo3(get)]
    weights: Vec<f64>,
    #[pyo3(get)]
    weight_std_errors: Vec<f64>,
    #[pyo3(get)]
    chi_sq: f64,
    #[pyo3(get)]
    df: usize,
    #[pyo3(get)]
    p_value: f64,
    #[pyo3(get)]
    constrained: bool,
    /// Fit of every subset of the sources (admixtools2 ``popdrop``).
    #[pyo3(get)]
    popdrop: Vec<PyPopDropRow>,
}

/// One row of qpAdm ``popdrop``: the model fitted to a subset of the sources.
///
/// ``pattern`` has one character per source, "1" if it is dropped; ``weights``
/// is None for dropped sources. ``best`` marks the nested path from the full
/// model, each step dropping one more source; ``dofdiff``, ``chisqdiff`` and
/// ``p_nested`` compare a step with the model it came from. With
/// ``compat="admixtools2-2.0.10"`` they reproduce admixtools2's columns
/// instead (DEVIATIONS.md, ADMIXTOOLS2-0001).
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPopDropRow {
    #[pyo3(get)]
    pattern: String,
    #[pyo3(get)]
    wt: usize,
    #[pyo3(get)]
    dof: usize,
    #[pyo3(get)]
    chi_sq: f64,
    #[pyo3(get)]
    p_value: f64,
    #[pyo3(get)]
    f4rank: usize,
    #[pyo3(get)]
    weights: Vec<Option<f64>>,
    #[pyo3(get)]
    feasible: bool,
    #[pyo3(get)]
    best: bool,
    #[pyo3(get)]
    dofdiff: Option<usize>,
    #[pyo3(get)]
    chisqdiff: Option<f64>,
    #[pyo3(get)]
    p_nested: Option<f64>,
}

impl From<popmix_core::QpAdmResult> for PyQpAdmResult {
    fn from(result: popmix_core::QpAdmResult) -> Self {
        PyQpAdmResult {
            target: result.target.to_string(),
            sources: result.sources.iter().map(|s| s.to_string()).collect(),
            weights: result.weights,
            weight_std_errors: result.weight_std_errors,
            chi_sq: result.chi_sq,
            df: result.df,
            p_value: result.p_value,
            constrained: result.constrained,
            popdrop: result
                .popdrop
                .into_iter()
                .map(|r| PyPopDropRow {
                    pattern: r.pattern,
                    wt: r.wt,
                    dof: r.dof,
                    chi_sq: r.chi_sq,
                    p_value: r.p_value,
                    f4rank: r.f4rank,
                    weights: r.weights,
                    feasible: r.feasible,
                    best: r.best,
                    dofdiff: r.dofdiff,
                    chisqdiff: r.chisqdiff,
                    p_nested: r.p_nested,
                })
                .collect(),
        }
    }
}

#[pymethods]
impl PyQpAdmResult {
    /// Check if the model is a good fit (p > threshold).
    #[pyo3(signature = (p_threshold=0.05))]
    fn is_good_fit(&self, p_threshold: f64) -> bool {
        self.p_value > p_threshold
    }

    /// Check if all weights are feasible (in [0, 1]).
    fn is_feasible(&self) -> bool {
        self.weights.iter().all(|&w| (0.0..=1.0).contains(&w))
    }

    /// Get weights as a dictionary mapping source name to weight, in source order.
    fn weights_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        for (source, weight) in self.sources.iter().zip(&self.weights) {
            dict.set_item(source, weight)?;
        }
        Ok(dict)
    }
}

/// A single model result from qpAdm rotation.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRotationModel {
    #[pyo3(get)]
    sources: Vec<String>,
    #[pyo3(get)]
    weights: Vec<f64>,
    #[pyo3(get)]
    weight_std_errors: Vec<f64>,
    #[pyo3(get)]
    p_value: f64,
    #[pyo3(get)]
    chi_sq: f64,
    #[pyo3(get)]
    df: usize,
    #[pyo3(get)]
    is_feasible: bool,
}

#[pymethods]
impl PyRotationModel {
    fn __repr__(&self) -> String {
        format!(
            "RotationModel(sources={:?}, p_value={:.4}, is_feasible={})",
            self.sources, self.p_value, self.is_feasible
        )
    }

    /// Get weights as a dictionary mapping source name to weight, in source order.
    fn weights_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        for (source, weight) in self.sources.iter().zip(&self.weights) {
            dict.set_item(source, weight)?;
        }
        Ok(dict)
    }
}

/// Result of qpAdm rotation analysis.
#[pyclass]
pub(crate) struct PyQpAdmRotateResult {
    #[pyo3(get)]
    target: String,
    #[pyo3(get)]
    candidates: Vec<String>,
    #[pyo3(get)]
    outgroups: Vec<String>,
    #[pyo3(get)]
    models: Vec<PyRotationModel>,
    #[pyo3(get)]
    n_combinations: usize,
    #[pyo3(get)]
    n_feasible: usize,
    /// Combinations not fit: more sources than outgroups, or a singular or
    /// infeasible system.
    #[pyo3(get)]
    n_skipped: usize,
}

#[pymethods]
impl PyQpAdmRotateResult {
    /// Get the best model (highest p-value).
    fn best_model(&self) -> Option<PyRotationModel> {
        self.models.first().cloned()
    }

    /// Get models passing a p-value threshold.
    #[pyo3(signature = (p_threshold=0.05))]
    fn passing_models(&self, p_threshold: f64) -> Vec<PyRotationModel> {
        self.models
            .iter()
            .filter(|m| m.p_value > p_threshold)
            .cloned()
            .collect()
    }

    /// Get feasible models only.
    fn feasible_models(&self) -> Vec<PyRotationModel> {
        self.models
            .iter()
            .filter(|m| m.is_feasible)
            .cloned()
            .collect()
    }
}

/// Run qpAdm rotation to test all source combinations.
///
/// Tests all combinations of candidate source populations to find
/// the best-fitting admixture models.
///
/// Args:
///     f2_blocks: Precomputed F2 blocks
///     target: Target population name
///     candidates: Candidate source population names
///     outgroups: Outgroup population names
///     min_sources: Minimum number of sources to test (default: 1)
///     max_sources: Maximum number of sources to test (default: all)
///     feasible_only: Only return feasible models (default: False)
///     fudge: Regularization factor (default: 1e-4)
///
/// Returns:
///     PyQpAdmRotateResult with all tested models sorted by p-value
#[pyfunction]
#[pyo3(signature = (f2_blocks, target, candidates, outgroups, min_sources=None, max_sources=None, feasible_only=None, fudge=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn qpadm_rotate(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    target: String,
    candidates: Vec<String>,
    outgroups: Vec<String>,
    min_sources: Option<usize>,
    max_sources: Option<usize>,
    feasible_only: Option<bool>,
    fudge: Option<f64>,
) -> PyResult<PyQpAdmRotateResult> {
    compute(py, move || {
        let candidates_refs: Vec<&str> = candidates.iter().map(|s| s.as_str()).collect();
        let outgroups_refs: Vec<&str> = outgroups.iter().map(|s| s.as_str()).collect();

        let opts = popmix_core::QpAdmRotateOptions {
            min_sources: min_sources.unwrap_or(1),
            max_sources,
            p_threshold: 0.05,
            feasible_only: feasible_only.unwrap_or(false),
            qpadm_opts: popmix_core::QpAdmOptions {
                fudge: fudge.unwrap_or(1e-4),
                ..Default::default()
            },
        };

        let result = popmix_core::qpadm_rotate(
            &f2_blocks.inner,
            &target,
            &candidates_refs,
            &outgroups_refs,
            Some(opts),
        )?;

        let models = result
            .models
            .into_iter()
            .map(|m| PyRotationModel {
                sources: m.sources.iter().map(|s| s.to_string()).collect(),
                weights: m.weights,
                weight_std_errors: m.weight_std_errors,
                p_value: m.p_value,
                chi_sq: m.chi_sq,
                df: m.df,
                is_feasible: m.is_feasible,
            })
            .collect();

        Ok(PyQpAdmRotateResult {
            target: result.target.to_string(),
            candidates: result.candidates.iter().map(|s| s.to_string()).collect(),
            outgroups: result.outgroups.iter().map(|s| s.to_string()).collect(),
            models,
            n_combinations: result.n_combinations,
            n_feasible: result.n_feasible,
            n_skipped: result.n_skipped,
        })
    })
}

/// Run qpWave analysis.
///
/// Tests whether left populations can be modeled as mixtures of ancestral
/// populations distinguished by the right (outgroup) populations.
///
/// Args:
///     f2_blocks: Precomputed F2 blocks
///     left: Left population names
///     right: Right (outgroup) population names
///     fudge: Regularization factor (default: 1e-4)
///
/// Returns:
///     PyQpWaveResult with rank tests and estimated rank
#[pyfunction]
#[pyo3(signature = (f2_blocks, left, right, fudge=None))]
pub(crate) fn qpwave(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    left: Vec<String>,
    right: Vec<String>,
    fudge: Option<f64>,
) -> PyResult<PyQpWaveResult> {
    compute(py, move || {
        let left_refs: Vec<&str> = left.iter().map(|s| s.as_str()).collect();
        let right_refs: Vec<&str> = right.iter().map(|s| s.as_str()).collect();

        let opts = popmix_core::QpWaveOptions {
            fudge: fudge.unwrap_or(1e-4),
            max_rank: None,
        };

        let result = popmix_core::qpwave(&f2_blocks.inner, &left_refs, &right_refs, Some(opts))?;

        let rank_tests = result
            .rank_tests
            .iter()
            .map(|r| PyRankTest {
                rank: r.rank,
                chi_sq: r.chi_sq,
                df: r.df,
                p_value: r.p_value,
            })
            .collect();

        let rows = |m: &ndarray::Array2<f64>| m.outer_iter().map(|r| r.to_vec()).collect();
        Ok(PyQpWaveResult {
            rank_tests,
            estimated_rank: result.estimated_rank,
            left_pops: left,
            right_pops: right,
            f4_matrix: rows(&result.f4_matrix.values),
            f4_covariance: rows(&result.f4_matrix.covariance),
        })
    })
}

/// Result of qpWave test for a single population pair.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQpWavePairResult {
    #[pyo3(get)]
    pop_a: String,
    #[pyo3(get)]
    pop_b: String,
    #[pyo3(get)]
    p_value: f64,
    #[pyo3(get)]
    chi_sq: f64,
    #[pyo3(get)]
    df: usize,
}

#[pymethods]
impl PyQpWavePairResult {
    fn __repr__(&self) -> String {
        format!(
            "QpWavePairResult(pop_a='{}', pop_b='{}', p_value={:.4})",
            self.pop_a, self.pop_b, self.p_value
        )
    }

    /// Check if this pair forms a clade (p > threshold).
    #[pyo3(signature = (p_threshold=0.05))]
    fn is_clade(&self, p_threshold: f64) -> bool {
        self.p_value > p_threshold
    }
}

impl From<&popmix_core::QpWavePairResult> for PyQpWavePairResult {
    fn from(p: &popmix_core::QpWavePairResult) -> Self {
        PyQpWavePairResult {
            pop_a: p.pop_a.to_string(),
            pop_b: p.pop_b.to_string(),
            p_value: p.p_value,
            chi_sq: p.chi_sq,
            df: p.df,
        }
    }
}

/// Result of pairwise qpWave analysis.
#[pyclass]
pub(crate) struct PyQpWavePairsResult {
    inner: popmix_core::QpWavePairsResult,
}

#[pymethods]
impl PyQpWavePairsResult {
    #[getter]
    fn populations(&self) -> Vec<String> {
        self.inner
            .populations
            .iter()
            .map(|s| s.to_string())
            .collect()
    }

    #[getter]
    fn outgroups(&self) -> Vec<String> {
        self.inner.outgroups.iter().map(|s| s.to_string()).collect()
    }

    #[getter]
    fn pairs(&self) -> Vec<PyQpWavePairResult> {
        self.inner.pairs.iter().map(Into::into).collect()
    }

    #[getter]
    fn n_pairs(&self) -> usize {
        self.inner.n_pairs
    }

    /// Pairs whose qpWave system was singular, left out of `pairs`.
    #[getter]
    fn n_skipped(&self) -> usize {
        self.inner.n_skipped
    }

    /// Get pairs that pass the cladality test (p > threshold).
    #[pyo3(signature = (p_threshold=0.05))]
    fn clades(&self, p_threshold: f64) -> Vec<PyQpWavePairResult> {
        self.inner
            .clades(p_threshold)
            .into_iter()
            .map(Into::into)
            .collect()
    }

    /// Get the pair with highest p-value (most likely clade); pairs with no
    /// p-value (0 degrees of freedom) are not candidates.
    fn best_clade(&self) -> Option<PyQpWavePairResult> {
        self.inner.best_clade().map(Into::into)
    }
}

/// Run qpWave for all population pairs (cladality screening).
///
/// Tests whether each pair of populations forms a clade relative to
/// the outgroup populations.
///
/// Args:
///     f2_blocks: Precomputed F2 blocks
///     populations: Population names to test pairwise
///     outgroups: Outgroup population names
///     fudge: Regularization factor (default: 1e-4)
///
/// Returns:
///     PyQpWavePairsResult with all pairwise test results
#[pyfunction]
#[pyo3(signature = (f2_blocks, populations, outgroups, fudge=None))]
pub(crate) fn qpwave_pairs(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    populations: Vec<String>,
    outgroups: Vec<String>,
    fudge: Option<f64>,
) -> PyResult<PyQpWavePairsResult> {
    compute(py, move || {
        let pop_refs: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
        let outgroup_refs: Vec<&str> = outgroups.iter().map(|s| s.as_str()).collect();

        let opts = popmix_core::QpWaveOptions {
            fudge: fudge.unwrap_or(1e-4),
            max_rank: None,
        };

        let inner =
            popmix_core::qpwave_pairs(&f2_blocks.inner, &pop_refs, &outgroup_refs, Some(opts))?;
        Ok(PyQpWavePairsResult { inner })
    })
}

/// Parse a qpAdm `compat` argument.
pub(crate) fn admixtools2_compat(compat: Option<&str>) -> PyResult<popmix_core::Admixtools2Compat> {
    match compat {
        None => Ok(popmix_core::Admixtools2Compat::CORRECTED),
        Some("admixtools2-2.0.10") => Ok(popmix_core::Admixtools2Compat::V2_0_10),
        Some(other) => Err(input_error(format!(
            "unknown compat {other:?} for qpadm; expected None or \"admixtools2-2.0.10\""
        ))),
    }
}

/// Run qpAdm analysis to estimate admixture weights.
///
/// Estimates the proportion of ancestry in a target population from
/// source populations, using outgroup populations to constrain the model.
///
/// Args:
///     f2_blocks: Precomputed F2 blocks
///     target: Target population name
///     sources: Source population names
///     right: Right (outgroup) population names
///     constrained: Whether to constrain weights to [0, 1] (default: False)
///     fudge: Regularization factor (default: 1e-4)
///     compat: Reproduce admixtools2's known quirks instead of correcting
///         them: "admixtools2-2.0.10", or None (default) for popmix's
///         corrected behaviour. See DEVIATIONS.md.
///
/// Returns:
///     PyQpAdmResult with estimated weights and model fit statistics
#[pyfunction]
#[pyo3(signature = (f2_blocks, target, sources, right, constrained=None, fudge=None, compat=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn qpadm(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    target: String,
    sources: Vec<String>,
    right: Vec<String>,
    constrained: Option<bool>,
    fudge: Option<f64>,
    compat: Option<&str>,
) -> PyResult<PyQpAdmResult> {
    compute(py, move || {
        let sources_refs: Vec<&str> = sources.iter().map(|s| s.as_str()).collect();
        let right_refs: Vec<&str> = right.iter().map(|s| s.as_str()).collect();

        let opts = popmix_core::QpAdmOptions {
            fudge: fudge.unwrap_or(1e-4),
            constrained: constrained.unwrap_or(false),
            compat: admixtools2_compat(compat)?,
        };

        let result = popmix_core::qpadm(
            &f2_blocks.inner,
            &target,
            &sources_refs,
            &right_refs,
            Some(opts),
        )?;

        Ok(PyQpAdmResult::from(result))
    })
}

/// Run qpAdm for multiple targets with the same sources and outgroups.
///
/// This is a convenience function for testing multiple target populations
/// against the same admixture model.
///
/// Args:
///     f2_blocks: Pre-computed F2 blocks
///     targets: List of target populations to test
///     sources: Source populations for the admixture model
///     right: Right (outgroup) populations
///     constrained: Use constrained optimization (weights in [0,1])
///     fudge: Regularization parameter
///     compat: "admixtools2-2.0.10" to reproduce admixtools2's quirks (see
///         qpadm)
///
/// Returns:
///     List of qpAdm results for each target
#[pyfunction]
#[pyo3(signature = (f2_blocks, targets, sources, right, constrained=None, fudge=None, compat=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn qpadm_multi(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    targets: Vec<String>,
    sources: Vec<String>,
    right: Vec<String>,
    constrained: Option<bool>,
    fudge: Option<f64>,
    compat: Option<&str>,
) -> PyResult<Vec<PyQpAdmResult>> {
    compute(py, move || {
        let sources_refs: Vec<&str> = sources.iter().map(|s| s.as_str()).collect();
        let right_refs: Vec<&str> = right.iter().map(|s| s.as_str()).collect();

        let opts = popmix_core::QpAdmOptions {
            fudge: fudge.unwrap_or(1e-4),
            constrained: constrained.unwrap_or(false),
            compat: admixtools2_compat(compat)?,
        };

        let mut results = Vec::with_capacity(targets.len());

        for target in &targets {
            match popmix_core::qpadm(
                &f2_blocks.inner,
                target,
                &sources_refs,
                &right_refs,
                Some(opts.clone()),
            ) {
                Ok(result) => {
                    results.push(PyQpAdmResult::from(result));
                }
                Err(e) => {
                    return Err(PyErr::from(e));
                }
            }
        }

        Ok(results)
    })
}

// ============================================================================
// qpGraph wrappers
// ============================================================================
