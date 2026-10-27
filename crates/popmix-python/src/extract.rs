//! `extract_f2`: f2 blocks computed from a genotype file, with population
//! relabelling and ploidy detection (`popmix_io::extract`).

use crate::compute;
use crate::f2::PyF2Blocks;
use crate::warnings::with_warnings;
use popmix_core::python::input_error;
use popmix_core::{BlockDefinition, F2Options};
use popmix_io::extract::{sample_labels, ExtractOptions};
use pyo3::prelude::*;

/// Extract F2 blocks from genotype data.
///
/// This is the main entry point for computing f-statistics.
/// Uses a streaming approach (like admixtools2) to avoid loading all genotype
/// data into memory. Reads SNPs one at a time and computes allele frequencies
/// incrementally.
///
/// Args:
///     prefix: Path prefix for the genotype dataset
///     populations: List of populations to compute F2 for (REQUIRED)
///     auto_only: If True (default), only use autosomal SNPs (chr 1-22)
///     block_size: Block size in Morgans (default: 0.05 = 5 cM)
///     maxmiss: Drop SNPs missing in more than this fraction of populations,
///         as admixtools2 ``maxmiss`` (default: 0, every population needs data)
///     poly_only: Use only SNPs whose frequency differs between populations,
///         as admixtools2 ``poly_only = "f2"`` (default: True)
///     inds, pops: Use only these samples, each in the population given in
///         ``pops`` (same length) instead of its ``.ind`` label, as
///         admixtools2 ``f2_from_geno(inds, pops)``
///     pseudohaploid: Which samples count one allele per SNP: "detect"
///         (default: those with no heterozygous call at any SNP),
///         "admixtools2" (no heterozygous call among the first 1,000 SNPs,
///         which also catches some diploid genomes; ADMIXTOOLS2-0007) or
///         "diploid" (none)
///
/// Returns:
///     PyF2Blocks object for computing f2, f3, f4 statistics
#[pyfunction]
#[pyo3(signature = (prefix, populations, auto_only=true, block_size=None, maxmiss=0.0, poly_only=true, inds=None, pops=None, pseudohaploid="detect"))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn extract_f2(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    auto_only: bool,
    block_size: Option<f64>,
    maxmiss: f64,
    poly_only: bool,
    inds: Option<Vec<String>>,
    pops: Option<Vec<String>>,
    pseudohaploid: &str,
) -> PyResult<PyF2Blocks> {
    let labels = match (inds, pops) {
        (None, None) => None,
        (Some(inds), Some(pops)) => Some(sample_labels(inds, pops)?),
        _ => return Err(input_error("give inds and pops together")),
    };
    let mut opts = ExtractOptions::new(populations);
    opts.auto_only = auto_only;
    if let Some(size) = block_size {
        opts.blocks = BlockDefinition::new(size);
    }
    opts.f2 = F2Options {
        max_missing: maxmiss,
        poly_only,
    };
    opts.labels = labels;
    opts.ploidy = pseudohaploid.parse()?;
    let inner = compute(py, || {
        with_warnings(|warn| popmix_io::extract::extract_f2(prefix, &opts, warn))
    })?;
    Ok(PyF2Blocks { inner })
}
