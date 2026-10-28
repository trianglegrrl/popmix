//! Python bindings for popmix.
//!
//! This crate provides PyO3-based Python bindings for the popmix library.

// PyO3 macro expansions generate these warnings
#![allow(clippy::useless_conversion)]

mod admixture;
mod convert;
mod coverage;
mod datasets;
mod dates;
mod extract;
mod f2;
mod geno_fstats;
mod genotypes;
mod kinship;
mod pca;
mod ploidy;
mod qpadm;
mod qpadm_geno;
mod qpadm_models;
mod qpgraph;
mod rcompat;
mod rolloff;
mod sfs;
mod subset;
mod warnings;

use popmix_core::python::popmix_exception;
use pyo3::prelude::*;

use crate::dates::{PyDatesResult, PyDecayFit, PyLdBin};
use crate::f2::PyF2Blocks;
use crate::genotypes::{PyAlleleFrequencies, PyGenotypeMatrix, PySample, PySnp};
use crate::qpadm::{PyF4RatioResult, PyPopDropRow, PyQpAdmResult, PyQpWaveResult, PyRankTest};
use crate::qpgraph::{
    PyAdmixGraph, PyQpGraphEdge, PyQpGraphF3, PyQpGraphResult, PyQpGraphStart, PyScoredGraph,
};
use convert::{convertf, merge};
use popmix_core::F2Blocks;

/// An error naming any of `pops` that is not in the f2 blocks.
fn require_populations(f2: &F2Blocks, pops: &[&str]) -> PyResult<()> {
    let unknown: Vec<&str> = pops
        .iter()
        .copied()
        .filter(|p| f2.pop_index(p).is_none())
        .collect();
    if unknown.is_empty() {
        return Ok(());
    }
    let have: Vec<String> = f2
        .populations()
        .iter()
        .map(std::string::ToString::to_string)
        .collect();
    Err(popmix_exception(
        "PopulationNotFoundError",
        format!(
            "not in the f2 blocks: '{}' (they have: {})",
            unknown.join("', '"),
            have.join(", ")
        ),
    ))
}

/// Requested populations each once, and each with samples in the data.
fn require_sample_populations(samples: &[popmix_core::Sample], pops: &[String]) -> PyResult<()> {
    let present: std::collections::HashSet<&str> =
        samples.iter().map(|s| s.population.as_str()).collect();
    popmix_core::genotypes::check_populations(pops.iter().map(std::string::String::as_str), |p| {
        present.contains(p)
    })
    .map_err(PyErr::from)
}

// ============================================================================
// Thread pool functions
// ============================================================================

/// Set the number of threads popmix computes with (4 until set).
///
/// It can be changed at any time.
///
/// Parameters
/// ----------
/// threads
///     Number of threads, from 1 to 1024; None (the default) means 4.
///
/// Returns
/// -------
/// bool
///     True. Raises ``InputError`` for a count outside 1 to 1024, and
///     ``FileError`` if the threads cannot be started.
#[pyfunction]
#[pyo3(signature = (threads=None))]
fn init_threads(threads: Option<usize>) -> PyResult<bool> {
    popmix_core::init_thread_pool(threads).map_err(PyErr::from)
}

/// Run `f` without the GIL on popmix's thread pool (`init_threads`).
pub(crate) fn compute<T: Send>(py: Python<'_>, f: impl FnOnce() -> T + Send) -> T {
    py.detach(|| popmix_core::threading::install(f))
}

/// Get the number of threads popmix computes with.
///
/// Returns
/// -------
/// int
///     The count set by ``init_threads``, or 4 if it has not been called.
#[pyfunction]
fn get_threads() -> usize {
    popmix_core::threading::current_threads()
}

/// Get the popmix version.
///
/// Returns
/// -------
/// str
///     The version of the compiled extension, such as ``"1.0.0"``.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

// ============================================================================
// Sample wrapper
// ============================================================================

// ============================================================================
// Module definition
// ============================================================================

/// Python module for popmix population genetics toolkit.
// Not yet tested on free-threaded Python (3.13t): until it is, importing
// popmix there turns the GIL back on rather than risk a race.
#[pymodule(gil_used = true)]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Threading functions
    // rebuilders for pickling
    m.add_function(wrap_pyfunction!(genotypes::_snp, m)?)?;
    m.add_function(wrap_pyfunction!(genotypes::_allele_frequencies, m)?)?;
    m.add_function(wrap_pyfunction!(f2::_f2_blocks, m)?)?;
    m.add_function(wrap_pyfunction!(init_threads, m)?)?;
    m.add_function(wrap_pyfunction!(get_threads, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;

    // I/O functions
    m.add_function(wrap_pyfunction!(genotypes::read_plink, m)?)?;
    m.add_function(wrap_pyfunction!(genotypes::read_eigenstrat, m)?)?;
    m.add_function(wrap_pyfunction!(genotypes::read_packedancestrymap, m)?)?;
    m.add_function(wrap_pyfunction!(genotypes::read_vcf, m)?)?;
    m.add_function(wrap_pyfunction!(genotypes::read_genotypes, m)?)?;
    m.add_function(wrap_pyfunction!(f2::read_f2, m)?)?;
    m.add_function(wrap_pyfunction!(datasets::dataset_info, m)?)?;
    m.add_function(wrap_pyfunction!(datasets::list_populations, m)?)?;
    m.add_function(wrap_pyfunction!(datasets::suggest_populations, m)?)?;
    m.add_function(wrap_pyfunction!(convertf, m)?)?;
    m.add_function(wrap_pyfunction!(merge, m)?)?;
    m.add_function(wrap_pyfunction!(datasets::detect_format, m)?)?;

    // F-statistics functions
    m.add_function(wrap_pyfunction!(extract::extract_f2, m)?)?;
    m.add_function(wrap_pyfunction!(f2::f2_from_frequencies, m)?)?;
    geno_fstats::register(m)?;
    pca::register(m)?;
    kinship::register(m)?;
    sfs::register(m)?;
    admixture::register(m)?;
    qpadm_models::register(m)?;
    rcompat::register(m)?;
    rolloff::register(m)?;
    subset::register(m)?;
    qpadm_geno::register(m)?;
    coverage::register(m)?;

    // qpAdm/qpWave functions
    m.add_function(wrap_pyfunction!(qpadm::qpwave, m)?)?;
    m.add_function(wrap_pyfunction!(qpadm::qpadm, m)?)?;
    m.add_function(wrap_pyfunction!(qpadm::qpf4ratio, m)?)?;

    // qpGraph functions
    m.add_function(wrap_pyfunction!(qpgraph::build_tree, m)?)?;
    m.add_function(wrap_pyfunction!(qpgraph::qpgraph, m)?)?;
    m.add_function(wrap_pyfunction!(qpgraph::find_graphs, m)?)?;
    m.add_function(wrap_pyfunction!(qpgraph::all_graphs, m)?)?;
    m.add_function(wrap_pyfunction!(qpgraph::graph_key, m)?)?;

    // DATES functions
    m.add_function(wrap_pyfunction!(dates::dates, m)?)?;

    // Classes
    m.add_class::<PySample>()?;
    m.add_class::<PySnp>()?;
    m.add_class::<PyGenotypeMatrix>()?;
    m.add_class::<PyAlleleFrequencies>()?;
    m.add_class::<PyF2Blocks>()?;
    m.add_class::<PyRankTest>()?;
    m.add_class::<PyQpWaveResult>()?;
    m.add_class::<PyQpAdmResult>()?;
    m.add_class::<PyPopDropRow>()?;
    m.add_class::<PyQpGraphResult>()?;
    m.add_class::<PyQpGraphEdge>()?;
    m.add_class::<PyQpGraphF3>()?;
    m.add_class::<PyQpGraphStart>()?;
    m.add_class::<PyF4RatioResult>()?;
    m.add_class::<PyAdmixGraph>()?;
    m.add_class::<PyScoredGraph>()?;
    m.add_class::<PyLdBin>()?;
    m.add_class::<PyDecayFit>()?;
    m.add_class::<PyDatesResult>()?;
    m.add_class::<dates::PyPulseCheck>()?;

    // Add version as module attribute
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
