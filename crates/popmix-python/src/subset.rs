//! Choosing and relabelling samples and SNPs (`popmix_io::subset`):
//! `convertf`'s selection, and the samples of f2 extraction on the CLI
//! (TASKS.md H2, P3).

use std::collections::{HashMap, HashSet};

use popmix_core::genotypes::Chromosome;
use popmix_io::subset::Subset;
use pyo3::prelude::*;

use crate::compute;
use crate::warnings::with_warnings;

/// A [`Subset`] from the keyword arguments of `convertf` and
/// `select_samples`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn subset(
    keep_inds: Option<Vec<String>>,
    exclude_inds: Option<Vec<String>>,
    relabel: Option<HashMap<String, String>>,
    relabel_pops: Option<HashMap<String, String>>,
    populations: Option<Vec<String>>,
    keep_snps: Option<Vec<String>>,
    exclude_snps: Option<Vec<String>>,
    chromosomes: Option<Vec<String>>,
    region: Option<(String, u32, u32)>,
) -> Subset {
    let set = |v: Option<Vec<String>>| v.map(|v| v.into_iter().collect::<HashSet<_>>());
    let mut s = Subset::default();
    s.keep_samples = set(keep_inds);
    s.exclude_samples = set(exclude_inds);
    s.relabel_samples = relabel.unwrap_or_default();
    s.relabel_populations = relabel_pops.unwrap_or_default();
    s.populations = populations;
    s.keep_snps = set(keep_snps);
    s.exclude_snps = set(exclude_snps);
    s.chromosomes = chromosomes.map(|c| c.iter().map(|x| Chromosome::parse(x)).collect());
    s.region = region.map(|(c, from, to)| (Chromosome::parse(&c), from, to));
    s
}

/// The samples of a dataset kept and relabelled, as ``extract_f2``'s
/// ``inds`` and ``pops``.
///
/// Samples are kept or excluded by ID, then relabelled by ID, then by
/// population. An ID or population not in the dataset raises
/// ``InputError``.
///
/// Parameters
/// ----------
/// prefix
///     Path prefix of the dataset, or the full path of a VCF file.
/// keep_inds, exclude_inds
///     Only these sample IDs; never these.
/// relabel
///     A sample ID's new population.
/// relabel_pops
///     A population's new name.
///
/// Returns
/// -------
/// tuple of list of str
///     The kept sample IDs, in file order, and their populations.
#[pyfunction]
#[pyo3(signature = (prefix, *, keep_inds=None, exclude_inds=None, relabel=None, relabel_pops=None))]
pub(crate) fn select_samples(
    py: Python<'_>,
    prefix: &str,
    keep_inds: Option<Vec<String>>,
    exclude_inds: Option<Vec<String>>,
    relabel: Option<HashMap<String, String>>,
    relabel_pops: Option<HashMap<String, String>>,
) -> PyResult<(Vec<String>, Vec<String>)> {
    let s = subset(
        keep_inds,
        exclude_inds,
        relabel,
        relabel_pops,
        None,
        None,
        None,
        None,
        None,
    );
    compute(py, || {
        let samples = with_warnings(|warn| popmix_io::dataset::read_samples(prefix, warn))?;
        let (_, kept) = s.samples(&samples)?;
        Ok(kept
            .into_iter()
            .map(|k| (k.id, k.population.as_str().to_string()))
            .unzip())
    })
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(select_samples, m)?)?;
    Ok(())
}
