//! Bindings for site-frequency spectra.

use ndarray::Array2;
use numpy::IntoPyArray;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use popmix_core::{joint_sfs, projected_sfs, AlleleFrequencies, JointSfsCompat, Spectrum};

use crate::ploidy::{ploidy_rule, PloidyRule};
use crate::to_py_err;

/// Allele frequencies and counts of `populations`, in that order, with
/// pseudo-haploid samples (found by `rule`) counting one allele.
fn allele_frequencies(
    prefix: &str,
    populations: &[String],
    auto_only: bool,
    rule: PloidyRule,
) -> PyResult<AlleleFrequencies> {
    let pops: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
    let af = if let Some(mut m) = crate::ploidy::mapped(prefix)? {
        rule.apply_mapped(&m.file, &mut m.samples, &pops);
        m.file
            .allele_frequencies(&m.samples, &m.snps, &pops, auto_only)
    } else {
        let mut geno = crate::read_genotypes_of(prefix, populations)?;
        rule.apply(&mut geno);
        let geno = if auto_only {
            geno.filter_autosomal().map_err(to_py_err)?
        } else {
            geno
        };
        AlleleFrequencies::from_genotypes_for_populations(&geno, Some(&pops))
    };
    // the requested order, and every population present
    let mut cols = Vec::with_capacity(pops.len());
    for p in &pops {
        let k = af
            .pop_index(p)
            .filter(|&k| af.counts.column(k).iter().any(|&c| c > 0));
        cols.push(k.ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "population {p} not in {prefix}"
            ))
        })?);
    }
    let n = af.frequencies.nrows();
    Ok(AlleleFrequencies {
        frequencies: Array2::from_shape_fn((n, cols.len()), |(s, j)| af.frequencies[[s, cols[j]]]),
        counts: Array2::from_shape_fn((n, cols.len()), |(s, j)| af.counts[[s, cols[j]]]),
        populations: cols.iter().map(|&k| af.populations[k].clone()).collect(),
        snps: af.snps,
    })
}

fn spectrum_dict<'py>(py: Python<'py>, sfs: Spectrum) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let dims = sfs.dims.clone();
    d.set_item("dims", dims.clone())?;
    let data = ndarray::ArrayD::from_shape_vec(dims.clone(), sfs.data)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
    let mask = ndarray::ArrayD::from_shape_vec(dims, sfs.mask)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
    d.set_item("data", data.into_pyarray(py))?;
    d.set_item("mask", mask.into_pyarray(py))?;
    Ok(d)
}

/// Joint SFS as admixtools2 ``joint_sfs()``: counts of SNPs by the allele
/// count of each population, from 0 to its largest sample size.
///
/// Returns a dict: dims, data (counts, first population slowest), mask, and
/// n_snps, in_sfs, missing, not_integer (what happened to the SNPs).
#[pyfunction(name = "joint_sfs")]
#[pyo3(signature = (prefix, populations, compat=None, auto_only=false))]
pub fn joint_sfs_py(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    compat: Option<&str>,
    auto_only: bool,
) -> PyResult<Py<PyAny>> {
    let compat = match compat {
        None => JointSfsCompat::CORRECTED,
        Some("admixtools2-2.0.10") => JointSfsCompat::ADMIXTOOLS2_2_0_10,
        Some(other) => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "unknown compat {other:?} for joint_sfs; expected None or \"admixtools2-2.0.10\""
            )))
        }
    };
    let af = crate::compute(py, || {
        allele_frequencies(prefix, &populations, auto_only, PloidyRule::First1000)
    })?;
    let (sfs, snps) = joint_sfs(&af, compat);
    let d = spectrum_dict(py, sfs)?;
    d.set_item("n_snps", snps.n_snps)?;
    d.set_item("in_sfs", snps.in_sfs)?;
    d.set_item("missing", snps.missing)?;
    d.set_item("not_integer", snps.not_integer)?;
    Ok(d.into())
}

/// SFS projected to ``sizes`` (allele counts per population), as dadi
/// ``Spectrum.from_data_dict``, folded if ``folded``.
///
/// Returns a dict: dims, data, mask (dadi's: corners, and folded-out entries)
/// and not_integer (SNPs left out for half an allele count).
#[pyfunction(name = "projected_sfs")]
#[pyo3(signature = (prefix, populations, sizes, folded=false, auto_only=true, pseudohaploid="detect"))]
pub fn projected_sfs_py(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    sizes: Vec<usize>,
    folded: bool,
    auto_only: bool,
    pseudohaploid: &str,
) -> PyResult<Py<PyAny>> {
    let rule = ploidy_rule(pseudohaploid)?;
    let af = crate::compute(py, || {
        allele_frequencies(prefix, &populations, auto_only, rule)
    })?;
    let max: Vec<usize> = (0..af.counts.ncols())
        .map(|p| af.counts.column(p).iter().copied().max().unwrap_or(0) as usize)
        .collect();
    let sizes = if sizes.is_empty() { max } else { sizes };
    let (sfs, not_integer) = projected_sfs(&af, &sizes, folded).map_err(to_py_err)?;
    let d = spectrum_dict(py, sfs)?;
    d.set_item("not_integer", not_integer)?;
    Ok(d.into())
}

/// dadi's statistics of a one-population spectrum: (S, pi, Watterson's
/// theta, Tajima's D).
#[pyfunction]
pub fn sfs_summary(data: Vec<f64>, mask: Vec<bool>) -> PyResult<(f64, f64, f64, f64)> {
    let sfs = Spectrum {
        dims: vec![data.len()],
        data,
        mask,
    };
    let s = popmix_core::SfsSummary::of(&sfs).map_err(to_py_err)?;
    Ok((s.s, s.pi, s.watterson_theta, s.tajima_d))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sfs_summary, m)?)?;
    m.add_function(wrap_pyfunction!(joint_sfs_py, m)?)?;
    m.add_function(wrap_pyfunction!(projected_sfs_py, m)?)?;
    Ok(())
}
