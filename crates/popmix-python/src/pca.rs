//! Bindings for smartpca-style PCA.

use numpy::{IntoPyArray, PyArray2};
use popmix_core::python::{input_error, internal_error};
use pyo3::prelude::*;

use popmix_core::{smartpca, EigensoftCompat, PackedGenotypes, PcaOptions, SampleRole};

use crate::geno_fstats::{open_source, rows_of};

/// (eigenvalue, twstat, p, neff) of one Tracy-Widom row.
type TwRow = (f64, Option<f64>, Option<f64>, Option<f64>);

/// Result of ``smartpca``: PCA samples and projected samples together.
#[pyclass(module = "popmix._core")]
pub struct PyPcaResult {
    /// Every eigenvalue, decreasing (smartpca's evaloutname).
    #[pyo3(get)]
    eigenvalues: Vec<f64>,
    /// Samples with coordinates, in .ind order: id, population, and
    /// whether it defines the PCA ("pca"), is projected ("project") or is a
    /// PCA population's sample held out and projected ("holdout").
    #[pyo3(get)]
    sample_ids: Vec<String>,
    #[pyo3(get)]
    populations: Vec<String>,
    #[pyo3(get)]
    roles: Vec<String>,
    coords: Vec<Vec<f64>>,
    /// Tracy-Widom table: (eigenvalue, twstat, p, neff), None where too few
    /// eigenvalues are left for a statistic.
    #[pyo3(get)]
    tracy_widom: Vec<TwRow>,
    /// Outliers removed: (id, population, iteration from 1, eigenvector from
    /// 0, score).
    #[pyo3(get)]
    outliers: Vec<(String, String, usize, usize, f64)>,
    /// Samples left out for lack of data.
    #[pyo3(get)]
    no_data: Vec<String>,
    #[pyo3(get)]
    n_pca_samples: usize,
    #[pyo3(get)]
    n_snps_used: usize,
    #[pyo3(get)]
    n_snps_killed: usize,
    /// Each sample's share of the PCA's SNPs it misses, as sample_ids.
    #[pyo3(get)]
    missing_fraction: Vec<f64>,
}

#[pymethods]
impl PyPcaResult {
    /// Coordinates, samples x components.
    #[getter]
    fn coords<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let k = self.coords.first().map_or(0, std::vec::Vec::len);
        let flat: Vec<f64> = self.coords.iter().flatten().copied().collect();
        let arr = ndarray::Array2::from_shape_vec((self.coords.len(), k), flat)
            .map_err(|e| internal_error(e.to_string()))?;
        Ok(arr.into_pyarray(py))
    }
}

/// PCA as EIGENSOFT smartpca: computed on the samples of ``pca_pops``, with
/// the samples of ``project_pops`` projected onto it.
///
/// Args:
///     prefix: Genotype prefix
///     pca_pops: Populations the PCA is computed on
///     project_pops: Populations projected onto it
///     n_components: Eigenvectors reported (smartpca numoutevec)
///     outlier_iterations: Rounds of outlier removal (numoutlieriter)
///     outlier_components: Eigenvectors examined for outliers (numoutlierevec)
///     outlier_sigma: Outlier threshold in standard deviations
///     lsqproject: Project by least squares on observed genotypes
///     min_tw_eigs: Fewest eigenvalues for a Tracy-Widom statistic (nostatslim)
///     holdout: Project every second sample of each PCA population
///     compat: None (corrected) or "eigensoft-7.2.1"
///     auto_only: Use only autosomes 1-22, as smartpca
#[pyfunction(name = "smartpca")]
#[pyo3(signature = (prefix, pca_pops, project_pops=vec![], n_components=10, outlier_iterations=5,
                    outlier_components=10, outlier_sigma=6.0, lsqproject=false, min_tw_eigs=10,
                    compat=None, auto_only=true, holdout=false))]
#[allow(clippy::too_many_arguments)]
pub fn smartpca_py(
    py: Python<'_>,
    prefix: &str,
    pca_pops: Vec<String>,
    project_pops: Vec<String>,
    n_components: usize,
    outlier_iterations: usize,
    outlier_components: usize,
    outlier_sigma: f64,
    lsqproject: bool,
    min_tw_eigs: usize,
    compat: Option<&str>,
    auto_only: bool,
    holdout: bool,
) -> PyResult<PyPcaResult> {
    let compat = match compat {
        None => EigensoftCompat::CORRECTED,
        Some("eigensoft-7.2.1") => EigensoftCompat::V7_2_1,
        Some(other) => {
            return Err(input_error(format!(
                "unknown compat {other:?} for smartpca; expected None or \"eigensoft-7.2.1\""
            )))
        }
    };
    if let Some(p) = pca_pops.iter().find(|p| project_pops.contains(p)) {
        return Err(input_error(format!(
            "{p} is both a PCA and a projected population"
        )));
    }
    let populations: Vec<String> = pca_pops.iter().chain(&project_pops).cloned().collect();
    let (mut source, snps, pops, samples) = open_source(prefix, &populations, auto_only)?;
    let held = held_out(&pops.sample_pop, pca_pops.len(), holdout)?;
    let roles: Vec<SampleRole> = pops
        .sample_pop
        .iter()
        .zip(&held)
        .map(|(&p, &h)| {
            if p < pca_pops.len() && !h {
                SampleRole::Pca
            } else {
                SampleRole::Project
            }
        })
        .collect();
    let options = {
        let mut o = PcaOptions::default();
        o.n_components = n_components;
        o.outlier_iterations = outlier_iterations;
        o.outlier_components = outlier_components;
        o.outlier_sigma = outlier_sigma;
        o.lsq_project = lsqproject;
        o.min_tw_eigs = min_tw_eigs;
        o.compat = compat;
        o
    };
    let n_samples = samples.len();
    let out = crate::compute(py, || {
        let mut rows = rows_of(&mut source);
        let geno = PackedGenotypes::read(rows.as_mut(), n_samples, snps.len())?;
        smartpca(&geno, &roles, &options)
    })?;

    let in_pca: Vec<f64> = out
        .samples
        .iter()
        .zip(&out.missing_fraction)
        .filter(|(&s, _)| roles[s] == SampleRole::Pca)
        .map(|(_, &f)| f)
        .collect();
    crate::warnings::emit(|w| popmix_core::warnings::check_pca_sample_missing(&in_pca, w))?;
    if !lsqproject {
        let projected: Vec<f64> = out
            .samples
            .iter()
            .zip(&out.missing_fraction)
            .filter(|(&s, _)| roles[s] == SampleRole::Project)
            .map(|(_, &f)| f)
            .collect();
        crate::warnings::emit(|w| popmix_core::warnings::check_projection_missing(&projected, w))?;
    }
    let id = |s: usize| samples[s].id.clone();
    let pop = |s: usize| samples[s].population.as_str().to_string();
    Ok(PyPcaResult {
        eigenvalues: out.eigenvalues,
        sample_ids: out.samples.iter().map(|&s| id(s)).collect(),
        populations: out.samples.iter().map(|&s| pop(s)).collect(),
        roles: out
            .samples
            .iter()
            .map(|&s| match roles[s] {
                SampleRole::Pca => "pca".to_string(),
                SampleRole::Project if held[s] => "holdout".to_string(),
                SampleRole::Project => "project".to_string(),
            })
            .collect(),
        coords: out.coords,
        tracy_widom: out
            .tracy_widom
            .iter()
            .map(|r| (r.eigenvalue, r.twstat, r.p, r.neff))
            .collect(),
        outliers: out
            .outliers
            .iter()
            .map(|o| {
                (
                    id(o.sample),
                    pop(o.sample),
                    o.iteration,
                    o.component,
                    o.score,
                )
            })
            .collect(),
        no_data: out.no_data.iter().map(|&s| id(s)).collect(),
        n_pca_samples: out.n_pca_samples,
        n_snps_used: out.n_snps_used,
        missing_fraction: out.missing_fraction,
        n_snps_killed: out.n_snps_killed,
    })
}

/// With `holdout`, which samples are held out of the PCA: the 2nd, 4th, ...
/// sample of each PCA population (`sample_pop` below `n_pca`), in .ind
/// order. Measuring their shrinkage needs two populations with a sample to
/// hold out.
fn held_out(sample_pop: &[usize], n_pca: usize, holdout: bool) -> PyResult<Vec<bool>> {
    if !holdout {
        return Ok(vec![false; sample_pop.len()]);
    }
    let mut seen = vec![0usize; n_pca];
    let held: Vec<bool> = sample_pop
        .iter()
        .map(|&p| {
            p < n_pca && {
                seen[p] += 1;
                seen[p].is_multiple_of(2)
            }
        })
        .collect();
    let split = seen.iter().filter(|&&n| n >= 2).count();
    if split < 2 {
        return Err(input_error(format!(
            "holdout needs two PCA populations with at least 2 samples each, to measure shrinkage; \
             {split} here"
        )));
    }
    Ok(held)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(smartpca_py, m)?)?;
    m.add_class::<PyPcaResult>()?;
    Ok(())
}
