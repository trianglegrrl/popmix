//! Bindings for ADMIXTURE-style ancestry estimation.

use numpy::{IntoPyArray, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use popmix_core::{AdmixtureOptions, PackedGenotypes, Sample, Snp};

use crate::geno_fstats::{open_source_filtered, rows_of};
use crate::to_py_err;

/// Genotypes of the samples of `populations`, with those samples and the
/// SNPs used: the SNPs in `snp_ids` (in file order), or else the SNPs on
/// `chromosomes` (every SNP if None) that some sample has called.
fn load(
    py: Python<'_>,
    prefix: &str,
    populations: &[String],
    chromosomes: &Option<Vec<String>>,
    snp_ids: &Option<Vec<String>>,
) -> PyResult<(PackedGenotypes, Vec<Snp>, Vec<Sample>)> {
    let chroms: Option<Vec<popmix_core::Chromosome>> = chromosomes.as_ref().map(|c| {
        c.iter()
            .map(|x| popmix_core::Chromosome::parse(x))
            .collect()
    });
    let wanted: Option<std::collections::HashSet<&str>> = snp_ids
        .as_ref()
        .map(|v| v.iter().map(|s| s.as_str()).collect());
    let keep = |s: &Snp| {
        chroms.as_ref().is_none_or(|c| c.contains(&s.chromosome))
            && wanted.as_ref().is_none_or(|w| w.contains(s.id.as_str()))
    };
    let (mut source, snps, _, samples) = open_source_filtered(prefix, populations, &keep)?;
    let n = samples.len();
    let geno = crate::compute(py, || {
        let mut rows = rows_of(&mut source);
        PackedGenotypes::read(rows.as_mut(), n, snps.len())
    })
    .map_err(to_py_err)?;
    if let Some(w) = &wanted {
        if snps.len() != w.len() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "{} of the {} SNPs given are in {prefix}",
                snps.len(),
                w.len()
            )));
        }
        return Ok((geno, snps, samples));
    }
    let called: Vec<usize> = (0..geno.n_snps()).filter(|&j| geno.any_called(j)).collect();
    let snps = called.iter().map(|&j| snps[j].clone()).collect();
    Ok((geno.select_snps(&called), snps, samples))
}

fn result(
    py: Python<'_>,
    fit: popmix_core::AdmixtureFit,
    snps: &[Snp],
    samples: &[Sample],
    cv_error: Option<f64>,
) -> PyResult<Py<PyAny>> {
    let d = PyDict::new(py);
    let k = fit.k;
    let as2d = |v: Vec<f64>, rows: usize| {
        ndarray::Array2::from_shape_vec((rows, k), v)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))
    };
    d.set_item(
        "ids",
        samples.iter().map(|s| s.id.clone()).collect::<Vec<_>>(),
    )?;
    d.set_item(
        "pops",
        samples
            .iter()
            .map(|s| s.population.as_str().to_string())
            .collect::<Vec<_>>(),
    )?;
    d.set_item(
        "snps",
        snps.iter().map(|s| s.id.clone()).collect::<Vec<_>>(),
    )?;
    d.set_item("Q", as2d(fit.q, samples.len())?.into_pyarray(py))?;
    d.set_item("P", as2d(fit.p, snps.len())?.into_pyarray(py))?;
    d.set_item("loglik", fit.loglik)?;
    d.set_item("iterations", fit.iterations)?;
    d.set_item("converged", fit.converged)?;
    d.set_item("cv_error", cv_error)?;
    Ok(d.into())
}

/// ADMIXTURE-style ancestry proportions with K populations.
///
/// Modes: unsupervised; ``supervised`` (reference populations whose samples
/// keep Q at their population, clusters numbered in order of first
/// appearance as ADMIXTURE); or projection with ``fixed_p`` (SNPs x K, the
/// P of an earlier fit on the same SNPs).
///
/// Returns a dict: ids, pops, snps, Q, P, loglik, iterations, converged,
/// cv_error (with ``cv`` folds).
#[pyfunction(name = "admixture")]
#[pyo3(signature = (prefix, k, populations, supervised=None, fixed_p=None, cv=None, chromosomes=None,
                    snps=None, seed=1, tolerance=1e-4, max_iterations=10000))]
#[allow(clippy::too_many_arguments)]
pub fn admixture_py(
    py: Python<'_>,
    prefix: &str,
    k: usize,
    populations: Vec<String>,
    supervised: Option<Vec<String>>,
    fixed_p: Option<PyReadonlyArray2<'_, f64>>,
    cv: Option<usize>,
    chromosomes: Option<Vec<String>>,
    snps: Option<Vec<String>>,
    seed: u64,
    tolerance: f64,
    max_iterations: usize,
) -> PyResult<Py<PyAny>> {
    let (geno, snps, samples) = load(py, prefix, &populations, &chromosomes, &snps)?;
    let options = AdmixtureOptions {
        tolerance,
        max_iterations,
        seed,
    };
    let p_fixed: Option<Vec<f64>> = fixed_p.map(|p| p.as_array().iter().copied().collect());
    let (fit, cv_error) = crate::compute(py, || -> popmix_core::Result<_> {
        let fit = if let Some(p) = &p_fixed {
            popmix_core::admixture_projection(&geno, p, k, &options)?
        } else if let Some(refs) = &supervised {
            let mut order: Vec<&str> = Vec::new();
            for s in &samples {
                let pop = s.population.as_str();
                if refs.iter().any(|r| r == pop) && !order.contains(&pop) {
                    order.push(pop);
                }
            }
            if order.len() != k {
                return Err(popmix_core::PopmixError::DimensionMismatch {
                    context: "supervised ADMIXTURE".to_string(),
                    expected: format!("{k} reference populations, one per K"),
                    actual: order.len().to_string(),
                });
            }
            let labels: Vec<Option<usize>> = samples
                .iter()
                .map(|s| order.iter().position(|o| *o == s.population.as_str()))
                .collect();
            popmix_core::admixture_supervised(&geno, k, &labels, &options)?
        } else {
            popmix_core::admixture(&geno, k, &options)?
        };
        let cv_error = match cv {
            Some(f) => Some(popmix_core::admixture_cv(&geno, &fit, f, &options)?),
            None => None,
        };
        Ok((fit, cv_error))
    })
    .map_err(to_py_err)?;
    result(py, fit, &snps, &samples, cv_error)
}

/// Log-likelihood of Q (samples x K) and P (SNPs x K), as ADMIXTURE.
#[pyfunction]
#[pyo3(signature = (prefix, populations, q, p, chromosomes=None, snps=None))]
pub fn admixture_loglik(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    q: PyReadonlyArray2<'_, f64>,
    p: PyReadonlyArray2<'_, f64>,
    chromosomes: Option<Vec<String>>,
    snps: Option<Vec<String>>,
) -> PyResult<f64> {
    let (geno, _, _) = load(py, prefix, &populations, &chromosomes, &snps)?;
    let k = q.as_array().ncols();
    let q: Vec<f64> = q.as_array().iter().copied().collect();
    let p: Vec<f64> = p.as_array().iter().copied().collect();
    popmix_core::admixture_loglik(&geno, &q, &p, k).map_err(to_py_err)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(admixture_py, m)?)?;
    m.add_function(wrap_pyfunction!(admixture_loglik, m)?)?;
    Ok(())
}
