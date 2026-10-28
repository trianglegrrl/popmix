//! qpWave and qpAdm on genotypes, as admixtools2's `qpwave(prefix)` and
//! `qpadm(prefix)`: the f4 matrix from allele frequencies at the model's
//! SNPs (`popmix_core::f4_matrix_from_genotypes`, TASKS.md H11).

use popmix_core::python::input_error;
use pyo3::prelude::*;

use crate::compute;
use crate::geno_fstats::{open_source_with, run};
use crate::ploidy::ploidy_rule;
use crate::qpadm::{admixtools2_compat, wave_result, PyQpAdmResult, PyQpWaveResult};

/// The f4 matrix of `left` against `right` from the genotypes at `prefix`,
/// with the number of blocks with SNPs and of SNPs (the block weights'
/// sum).
fn f4_matrix(
    prefix: &str,
    left: &[String],
    right: &[String],
    allsnps: bool,
    block_size: f64,
    auto_only: bool,
    pseudohaploid: &str,
) -> PyResult<(popmix_core::F4Matrix, usize, usize)> {
    let mut populations: Vec<String> = Vec::new();
    for p in left.iter().chain(right) {
        if populations.contains(p) {
            return Err(input_error(format!("{p} is given twice")));
        }
        populations.push(p.clone());
    }
    let rule = ploidy_rule(pseudohaploid)?;
    let (source, snps, pops, _) = open_source_with(prefix, &populations, auto_only, rule)?;
    let index = |p: &String| {
        pops.names
            .iter()
            .position(|n| n == p)
            .expect("opened by name")
    };
    let (l, r): (Vec<usize>, Vec<usize>) = (
        left.iter().map(index).collect(),
        right.iter().map(index).collect(),
    );
    let (matrix, weights) = run(source, |rows| {
        popmix_core::f4_matrix_from_genotypes(rows, &snps, &pops, &l, &r, block_size, allsnps)
    })?;
    let n_blocks = weights.iter().filter(|&&w| w > 0.0).count();
    let statistics = (left.len() - 1) * (right.len() - 1);
    crate::warnings::emit(|w| popmix_core::warnings::check_blocks_for_f4(n_blocks, statistics, w))?;
    Ok((matrix, n_blocks, weights.iter().sum::<f64>() as usize))
}

/// qpAdm on genotypes (admixtools2 ``qpadm(prefix)``).
///
/// Args:
///     prefix: Genotype prefix
///     left, right, target, constrained, fudge, compat: as ``qpadm``
///     allsnps: Each f4 on the SNPs its four populations have data for,
///         rather than the SNPs every population of the model has
///     block_size: Jackknife block size in Morgans
///     auto_only: Use only autosomes 1-22
///     pseudohaploid: Ploidy rule ("detect", "admixtools2", "diploid")
///
/// Returns:
///     (PyQpAdmResult, blocks with SNPs, SNPs)
#[pyfunction]
#[pyo3(signature = (prefix, left, right, target, *, constrained=false, fudge=1e-4, compat=None,
                    allsnps=false, block_size=0.05, auto_only=true, pseudohaploid="detect"))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn qpadm_geno(
    py: Python<'_>,
    prefix: &str,
    left: Vec<String>,
    right: Vec<String>,
    target: String,
    constrained: bool,
    fudge: f64,
    compat: Option<&str>,
    allsnps: bool,
    block_size: f64,
    auto_only: bool,
    pseudohaploid: &str,
) -> PyResult<(PyQpAdmResult, usize, usize)> {
    let opts = {
        let mut o = popmix_core::QpAdmOptions::default();
        o.fudge = fudge;
        o.constrained = constrained;
        o.compat = admixtools2_compat(compat)?;
        o
    };
    compute(py, move || -> PyResult<_> {
        let mut wave_left = vec![target.clone()];
        wave_left.extend(left.iter().cloned());
        let (m, n_blocks, n_snps) = f4_matrix(
            prefix,
            &wave_left,
            &right,
            allsnps,
            block_size,
            auto_only,
            pseudohaploid,
        )?;
        let sources: Vec<&str> = left.iter().map(String::as_str).collect();
        let result = popmix_core::qpadm_on(&target, &sources, m, &opts)?;
        Ok((PyQpAdmResult::from(result), n_blocks, n_snps))
    })
}

/// qpWave on genotypes (admixtools2 ``qpwave(prefix)``).
///
/// Args:
///     prefix: Genotype prefix
///     left, right, fudge, alpha: as ``qpwave``
///     allsnps, block_size, auto_only, pseudohaploid: as ``qpadm_geno``
#[pyfunction]
#[pyo3(signature = (prefix, left, right, *, fudge=1e-4, alpha=0.05, allsnps=false, block_size=0.05,
                    auto_only=true, pseudohaploid="detect"))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn qpwave_geno(
    py: Python<'_>,
    prefix: &str,
    left: Vec<String>,
    right: Vec<String>,
    fudge: f64,
    alpha: f64,
    allsnps: bool,
    block_size: f64,
    auto_only: bool,
    pseudohaploid: &str,
) -> PyResult<PyQpWaveResult> {
    compute(py, move || -> PyResult<_> {
        let (m, _, _) = f4_matrix(
            prefix,
            &left,
            &right,
            allsnps,
            block_size,
            auto_only,
            pseudohaploid,
        )?;
        let opts = {
            let mut o = popmix_core::QpWaveOptions::default();
            o.fudge = fudge;
            o.alpha = alpha;
            o
        };
        let result = popmix_core::qpwave_on(m, &opts)?;
        Ok(wave_result(&result, alpha, left, right))
    })
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(qpadm_geno, m)?)?;
    m.add_function(wrap_pyfunction!(qpwave_geno, m)?)?;
    Ok(())
}
