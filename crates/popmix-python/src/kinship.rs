//! Bindings for kinship estimators.

use numpy::IntoPyArray;
use popmix_core::python::input_error;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use popmix_core::{king_robust, read_kinship, BitGenotypes, ReadCompat, ReadNorm, ReadOptions};

use crate::geno_fstats::{open_source, rows_of};

/// KING-robust kinship for every pair of samples of ``populations``, as
/// PLINK 2 ``--make-king-table counts``.
///
/// Returns a dict of arrays, one entry per pair in PLINK 2's order: id1,
/// id2 (sample IDs), nsnp, hethet, ibs0, het1_hom2, het2_hom1 (counts over
/// SNPs called in both) and kinship.
#[pyfunction]
#[pyo3(signature = (prefix, populations, auto_only=true))]
pub fn king(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    auto_only: bool,
) -> PyResult<Py<PyAny>> {
    let (mut source, snps, _, samples) = open_source(prefix, &populations, auto_only)?;
    let n = samples.len();
    let pairs = crate::compute(py, || {
        let mut rows = rows_of(&mut source);
        let g = BitGenotypes::read(rows.as_mut(), n, snps.len())?;
        Ok::<_, popmix_core::PopmixError>(king_robust(&g))
    })?;
    let d = PyDict::new(py);
    let id = |k: usize| samples[k].id.clone();
    d.set_item("id1", pairs.iter().map(|p| id(p.i)).collect::<Vec<_>>())?;
    d.set_item("id2", pairs.iter().map(|p| id(p.j)).collect::<Vec<_>>())?;
    let count = |f: fn(&popmix_core::KingPair) -> u64| pairs.iter().map(f).collect::<Vec<u64>>();
    d.set_item("nsnp", count(|p| p.nsnp).into_pyarray(py))?;
    d.set_item("hethet", count(|p| p.hethet).into_pyarray(py))?;
    d.set_item("ibs0", count(|p| p.ibs0).into_pyarray(py))?;
    d.set_item("het1_hom2", count(|p| p.het1_hom2).into_pyarray(py))?;
    d.set_item("het2_hom1", count(|p| p.het2_hom1).into_pyarray(py))?;
    d.set_item(
        "kinship",
        pairs
            .iter()
            .map(|p| p.kinship)
            .collect::<Vec<f64>>()
            .into_pyarray(py),
    )?;
    Ok(d.into())
}

/// READ-style kinship (READv2) for every pair of samples of ``populations``.
///
/// Returns a dict: norm_value, and arrays one entry per pair (``(0, 1), (0,
/// 2), ...`` in sample order, pairs with no SNP called in both left out):
/// id1, id2, relationship, z_upper, z_lower (NaN where READv2 has NA),
/// normalized_p0, normalized_p0_se, p0, p0_se, first_degree_type ("" for
/// none), first_degree_windows_outside, overlap_snps, effective_snps,
/// kinship, kinship_ci_low, kinship_ci_high.
#[pyfunction(name = "read_kinship")]
#[pyo3(signature = (prefix, populations, norm="median", norm_value=None, window_size=5_000_000,
                    window_estimate=false, thresholds="read", compat=None, auto_only=false))]
#[allow(clippy::too_many_arguments)]
pub fn read_kinship_py(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    norm: &str,
    norm_value: Option<f64>,
    window_size: u64,
    window_estimate: bool,
    thresholds: &str,
    compat: Option<&str>,
    auto_only: bool,
) -> PyResult<Py<PyAny>> {
    let bad = |msg: String| input_error(msg);
    let norm = match (norm, norm_value) {
        (_, Some(v)) => ReadNorm::Value(v),
        ("median", None) => ReadNorm::Median,
        ("mean", None) => ReadNorm::Mean,
        ("max", None) => ReadNorm::Max,
        (other, None) => {
            return Err(bad(format!(
                "norm must be median, mean or max (or give norm_value), got {other:?}"
            )))
        }
    };
    let thresholds = match thresholds {
        "read" => popmix_core::kinship::READ_THRESHOLDS,
        "2pow" => popmix_core::kinship::read_thresholds_2pow(),
        other => {
            return Err(bad(format!(
                "thresholds must be \"read\" or \"2pow\", got {other:?}"
            )))
        }
    };
    let compat = match compat {
        None => ReadCompat::CORRECTED,
        Some("readv2-2.01") => ReadCompat::V2_01,
        Some(other) => {
            return Err(bad(format!(
                "unknown compat {other:?} for READ; expected None or \"readv2-2.01\""
            )))
        }
    };
    let options = {
        let mut o = ReadOptions::default();
        o.window_size = window_size;
        o.norm = norm;
        o.window_estimate = window_estimate;
        o.thresholds = thresholds;
        o.compat = compat;
        o
    };
    let (mut source, snps, _, samples) = open_source(prefix, &populations, auto_only)?;
    let positions: Vec<u64> = snps.iter().map(|s| s.position_bp as u64).collect();
    let mut chroms = Vec::with_capacity(snps.len());
    let mut run = 0u32;
    for (k, s) in snps.iter().enumerate() {
        if k > 0 && s.chromosome != snps[k - 1].chromosome {
            run += 1;
        }
        chroms.push(run);
    }
    let n = samples.len();
    let out = crate::compute(py, || {
        let mut rows = rows_of(&mut source);
        let g = BitGenotypes::read(rows.as_mut(), n, snps.len())?;
        read_kinship(&g, &positions, &chroms, &options)
    })?;
    let d = PyDict::new(py);
    let p = &out.pairs;
    let id = |k: usize| samples[k].id.clone();
    let f =
        |g: &dyn Fn(&popmix_core::kinship::ReadPair) -> f64| p.iter().map(g).collect::<Vec<f64>>();
    d.set_item("norm_value", out.norm_value)?;
    d.set_item("id1", p.iter().map(|x| id(x.i)).collect::<Vec<_>>())?;
    d.set_item("id2", p.iter().map(|x| id(x.j)).collect::<Vec<_>>())?;
    d.set_item(
        "relationship",
        p.iter().map(|x| x.relationship.label()).collect::<Vec<_>>(),
    )?;
    d.set_item(
        "z_upper",
        f(&|x| x.z_upper.unwrap_or(f64::NAN)).into_pyarray(py),
    )?;
    d.set_item(
        "z_lower",
        f(&|x| x.z_lower.unwrap_or(f64::NAN)).into_pyarray(py),
    )?;
    d.set_item("normalized_p0", f(&|x| x.normalized_p0).into_pyarray(py))?;
    d.set_item(
        "normalized_p0_se",
        f(&|x| x.normalized_p0_se).into_pyarray(py),
    )?;
    d.set_item("p0", f(&|x| x.p0).into_pyarray(py))?;
    d.set_item("p0_se", f(&|x| x.p0_se).into_pyarray(py))?;
    d.set_item(
        "first_degree_type",
        p.iter()
            .map(|x| x.first_degree_type.map_or("", |t| t.label()))
            .collect::<Vec<_>>(),
    )?;
    d.set_item(
        "first_degree_windows_outside",
        f(&|x| x.first_degree_windows_outside).into_pyarray(py),
    )?;
    d.set_item(
        "overlap_snps",
        p.iter()
            .map(|x| x.overlap_snps)
            .collect::<Vec<u64>>()
            .into_pyarray(py),
    )?;
    d.set_item("effective_snps", f(&|x| x.effective_snps).into_pyarray(py))?;
    d.set_item("kinship", f(&|x| x.kinship).into_pyarray(py))?;
    d.set_item("kinship_ci_low", f(&|x| x.kinship_ci.0).into_pyarray(py))?;
    d.set_item("kinship_ci_high", f(&|x| x.kinship_ci.1).into_pyarray(py))?;
    Ok(d.into())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(king, m)?)?;
    m.add_function(wrap_pyfunction!(read_kinship_py, m)?)?;
    Ok(())
}
