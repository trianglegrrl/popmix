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
/// This is the main entry point for computing f-statistics. It streams the
/// SNPs one at a time (like admixtools2), computing allele frequencies
/// incrementally rather than loading all genotypes into memory.
///
/// Parameters
/// ----------
/// prefix
///     Path prefix of the genotype dataset.
/// populations
///     Populations to compute f2 for (required).
/// auto_only
///     If True (default), use only autosomal SNPs (chromosomes 1-22).
/// block_size
///     Jackknife block size in Morgans (default 0.05, 5 cM).
/// maxmiss
///     Drop SNPs missing in more than this fraction of populations, as
///     admixtools2 ``maxmiss`` (default 0: every population needs data).
/// poly_only
///     Use only SNPs whose frequency differs between populations, as
///     admixtools2 ``poly_only = "f2"`` (default True).
/// inds, pops
///     Use only these samples, each in the population given in ``pops``
///     (same length) instead of its ``.ind`` label, as admixtools2
///     ``f2_from_geno(inds, pops)``. Give both or neither (``InputError``).
/// pseudohaploid
///     Which samples count one allele per SNP: "detect" (default: those with
///     no heterozygous call at any SNP), "admixtools2" (no heterozygous call
///     among the first 1,000 SNPs, which also catches some diploid genomes;
///     ADMIXTOOLS2-0007) or "diploid" (none).
/// transitions, transversions
///     Keep transitions (A<->G, C<->T) and transversions (every other pair
///     of bases), as admixtools2; ``transitions=False`` for transversions
///     only, as ancient DNA without UDG treatment needs. A SNP without two
///     known bases is dropped by either.
/// minmaf, maxmaf
///     Drop SNPs whose minor allele frequency, over the populations'
///     alleles pooled, is below ``minmaf`` or above ``maxmaf`` (defaults 0
///     and 0.5), as admixtools2.
/// keepsnps
///     Only these SNP IDs, overriding ``maxmiss``, ``auto_only``, the MAF
///     bounds and the transition choice, as admixtools2 ``keepsnps``;
///     ``poly_only`` and ``excludesnps`` still apply.
/// excludesnps
///     Never these SNP IDs (popmix's own; admixtools2 has none).
///
/// Returns
/// -------
/// F2Blocks
///     The per-block f2 statistics of the populations, for f2, f3, f4 and
///     FST. Problems met while reading are issued as ``PopmixWarning``.
#[pyfunction]
#[pyo3(signature = (prefix, populations, *, auto_only=true, block_size=0.05, maxmiss=0.0, poly_only=true, inds=None, pops=None, pseudohaploid="detect", transitions=true, transversions=true, minmaf=0.0, maxmaf=0.5, keepsnps=None, excludesnps=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn extract_f2(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    auto_only: bool,
    block_size: f64,
    maxmiss: f64,
    poly_only: bool,
    inds: Option<Vec<String>>,
    pops: Option<Vec<String>>,
    pseudohaploid: &str,
    transitions: bool,
    transversions: bool,
    minmaf: f64,
    maxmaf: f64,
    keepsnps: Option<Vec<String>>,
    excludesnps: Option<Vec<String>>,
) -> PyResult<PyF2Blocks> {
    let labels = match (inds, pops) {
        (None, None) => None,
        (Some(inds), Some(pops)) => Some(sample_labels(inds, pops)?),
        _ => return Err(input_error("give inds and pops together")),
    };
    let mut opts = ExtractOptions::new(populations);
    opts.auto_only = auto_only;
    opts.blocks = BlockDefinition::new(block_size);
    opts.f2 = {
        let mut o = F2Options::default();
        o.max_missing = maxmiss;
        o.poly_only = poly_only;
        o.transitions = transitions;
        o.transversions = transversions;
        o.min_maf = minmaf;
        o.max_maf = maxmaf;
        o.keep_snps = keepsnps.map(|v| v.into_iter().collect());
        o.exclude_snps = excludesnps.map(|v| v.into_iter().collect());
        o
    };
    opts.labels = labels;
    opts.ploidy = pseudohaploid.parse()?;
    let inner = compute(py, || {
        with_warnings(|warn| popmix_io::extract::extract_f2(prefix, &opts, warn))
    })?;
    crate::warnings::checked_blocks(inner)
}
