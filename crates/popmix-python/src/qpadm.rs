//! qpWave, qpAdm and f4-ratio for Python: single models, rotations,
//! pairwise cladality and model sets.

use crate::compute;
use crate::f2::PyF2Blocks;
use popmix_core::python::input_error;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Result of a single rank test from qpWave.
#[pyclass(module = "popmix._core", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRankTest {
    #[pyo3(get)]
    rank: usize,
    #[pyo3(get, name = "chisq")]
    chi_sq: f64,
    #[pyo3(get, name = "dof")]
    df: usize,
    #[pyo3(get, name = "p")]
    p_value: f64,
}

/// Result of qpWave analysis.
#[pyclass(module = "popmix._core")]
pub(crate) struct PyQpWaveResult {
    #[pyo3(get)]
    rank_tests: Vec<PyRankTest>,
    #[pyo3(get)]
    estimated_rank: usize,
    /// The significance level of the rank tests.
    #[pyo3(get)]
    alpha: f64,
    #[pyo3(get, name = "left")]
    left_pops: Vec<String>,
    #[pyo3(get, name = "right")]
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
    /// Get the full rank test.
    fn full_rank(&self) -> Option<PyRankTest> {
        self.rank_tests.last().cloned()
    }
}

/// Result of f4 ratio analysis for admixture proportion estimation.
#[pyclass(module = "popmix._core")]
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
    #[pyo3(get, name = "se")]
    std_err: f64,
    #[pyo3(get, name = "z")]
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
            &Default::default(),
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
#[pyclass(module = "popmix._core")]
pub(crate) struct PyQpAdmResult {
    #[pyo3(get)]
    target: String,
    #[pyo3(get)]
    sources: Vec<String>,
    #[pyo3(get)]
    weights: Vec<f64>,
    #[pyo3(get, name = "weight_se")]
    weight_std_errors: Vec<f64>,
    #[pyo3(get, name = "chisq")]
    chi_sq: f64,
    #[pyo3(get, name = "dof")]
    df: usize,
    #[pyo3(get, name = "p")]
    p_value: f64,
    #[pyo3(get)]
    constrained: bool,
    /// Fit of every subset of the sources (admixtools2 ``popdrop``).
    #[pyo3(get)]
    popdrop: Vec<PyPopDropRow>,
    /// Rank tests of the f4 matrix of the target and sources against the
    /// right populations (admixtools2 ``rankdrop``).
    #[pyo3(get)]
    rankdrop: Vec<PyRankTest>,
}

/// One row of qpAdm ``popdrop``: the model fitted to a subset of the sources.
///
/// ``pattern`` has one character per source, "1" if it is dropped; ``weights``
/// is None for dropped sources. ``best`` marks the nested path from the full
/// model, each step dropping one more source; ``dofdiff``, ``chisqdiff`` and
/// ``p_nested`` compare a step with the model it came from. With
/// ``compat="admixtools2-2.0.10"`` they reproduce admixtools2's columns
/// instead (DEVIATIONS.md, ADMIXTOOLS2-0001).
#[pyclass(module = "popmix._core", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPopDropRow {
    #[pyo3(get)]
    pattern: String,
    #[pyo3(get)]
    wt: usize,
    #[pyo3(get)]
    dof: usize,
    #[pyo3(get, name = "chisq")]
    chi_sq: f64,
    #[pyo3(get, name = "p")]
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
            sources: result
                .sources
                .iter()
                .map(std::string::ToString::to_string)
                .collect(),
            weights: result.weights,
            weight_std_errors: result.weight_std_errors,
            chi_sq: result.chi_sq,
            df: result.df,
            p_value: result.p_value,
            constrained: result.constrained,
            rankdrop: result
                .qpwave_result
                .rank_tests
                .iter()
                .map(|r| PyRankTest {
                    rank: r.rank,
                    chi_sq: r.chi_sq,
                    df: r.df,
                    p_value: r.p_value,
                })
                .collect(),
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
    /// Whether every weight is feasible (in [0, 1]).
    #[getter]
    fn feasible(&self) -> bool {
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
#[pyo3(signature = (f2_blocks, left, right, *, fudge=1e-4, alpha=0.05))]
pub(crate) fn qpwave(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    left: Vec<String>,
    right: Vec<String>,
    fudge: f64,
    alpha: f64,
) -> PyResult<PyQpWaveResult> {
    let statistics = left.len().saturating_sub(1) * right.len().saturating_sub(1);
    let blocks = f2_blocks.inner.n_blocks();
    let out = compute(py, move || -> PyResult<_> {
        let left_refs: Vec<&str> = left.iter().map(std::string::String::as_str).collect();
        let right_refs: Vec<&str> = right.iter().map(std::string::String::as_str).collect();

        let opts = {
            let mut o = popmix_core::QpWaveOptions::default();
            o.fudge = fudge;
            o.alpha = alpha;
            o
        };

        let result = popmix_core::qpwave(&f2_blocks.inner, &left_refs, &right_refs, &opts)?;
        Ok(wave_result(&result, alpha, left, right))
    })?;
    crate::warnings::emit(|w| popmix_core::warnings::check_blocks_for_f4(blocks, statistics, w))?;
    Ok(out)
}

/// The Python qpWave result of `result`.
pub(crate) fn wave_result(
    result: &popmix_core::QpWaveResult,
    alpha: f64,
    left: Vec<String>,
    right: Vec<String>,
) -> PyQpWaveResult {
    {
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
        PyQpWaveResult {
            rank_tests,
            estimated_rank: result.estimated_rank,
            alpha,
            left_pops: left,
            right_pops: right,
            f4_matrix: rows(&result.f4_matrix.values),
            f4_covariance: rows(&result.f4_matrix.covariance),
        }
    }
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
///     left: Source population names (admixtools2's ``left``, the target
///         apart)
///     right: Right (outgroup) population names
///     target: Target population name
///     constrained: Whether to constrain weights to [0, 1] (default: False)
///     fudge: Regularization factor (default: 1e-4)
///     compat: Reproduce admixtools2's known quirks instead of correcting
///         them: "admixtools2-2.0.10", or None (default) for popmix's
///         corrected behaviour. See DEVIATIONS.md.
///
/// Returns:
///     PyQpAdmResult with estimated weights and model fit statistics
#[pyfunction]
#[pyo3(signature = (f2_blocks, left, right, target, *, constrained=false, fudge=1e-4, compat=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn qpadm(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    left: Vec<String>,
    right: Vec<String>,
    target: String,
    constrained: bool,
    fudge: f64,
    compat: Option<&str>,
) -> PyResult<PyQpAdmResult> {
    // f4(target, source; right0, right): (left - 1)(right - 1) statistics
    let statistics = left.len() * right.len().saturating_sub(1);
    let blocks = f2_blocks.inner.n_blocks();
    let out = compute(py, move || -> PyResult<_> {
        let sources_refs: Vec<&str> = left.iter().map(std::string::String::as_str).collect();
        let right_refs: Vec<&str> = right.iter().map(std::string::String::as_str).collect();

        let opts = {
            let mut o = popmix_core::QpAdmOptions::default();
            o.fudge = fudge;
            o.constrained = constrained;
            o.compat = admixtools2_compat(compat)?;
            o
        };

        let result =
            popmix_core::qpadm(&f2_blocks.inner, &target, &sources_refs, &right_refs, &opts)?;

        Ok(PyQpAdmResult::from(result))
    })?;
    crate::warnings::emit(|w| popmix_core::warnings::check_blocks_for_f4(blocks, statistics, w))?;
    Ok(out)
}
