//! Each sample's calls (`popmix_core::coverage`, TASKS.md H10).

use pyo3::prelude::*;

use crate::geno_fstats::{open_source, run};

/// (id, population, SNPs called, heterozygous calls).
type CallRow = (String, String, u64, u64);

/// Each sample's calls: (id, population, SNPs called, heterozygous calls),
/// in file order, and the number of SNPs read.
///
/// Args:
///     prefix: Genotype prefix
///     populations: Populations whose samples are counted
///     auto_only: Count only autosomal SNPs
#[pyfunction]
#[pyo3(signature = (prefix, populations, auto_only=true))]
pub(crate) fn sample_calls(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    auto_only: bool,
) -> PyResult<(Vec<CallRow>, usize)> {
    let (source, snps, _, samples) = open_source(prefix, &populations, auto_only)?;
    let (n, m) = (samples.len(), snps.len());
    let calls = crate::compute(py, || {
        run(source, |rows| {
            popmix_core::coverage::sample_calls(rows, n, m)
        })
    })?;
    let rows = samples
        .iter()
        .zip(calls)
        .map(|(s, c)| {
            (
                s.id.clone(),
                s.population.as_str().to_string(),
                c.called,
                c.heterozygous,
            )
        })
        .collect();
    Ok((rows, m))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sample_calls, m)?)?;
    Ok(())
}
