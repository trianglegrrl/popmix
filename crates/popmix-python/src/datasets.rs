//! What a dataset holds: its format, its size, its populations, and the
//! names nearest a mistyped one.

use crate::compute;
use crate::warnings::open_vcf;
use popmix_core::python::missing_file;
use popmix_core::Sample;
use popmix_io::dataset::vcf_path;
use popmix_io::tgeno::{is_tgeno_format, TgenoReader};
use popmix_io::{
    is_packed_format, EigenstratReader, GenoFormat, PackedAncestryMapReader, PlinkReader,
};
use pyo3::prelude::*;
use pyo3::IntoPyObjectExt;
use std::collections::HashMap;

/// Detect the genotype file format at a given prefix.
///
/// Parameters
/// ----------
/// prefix
///     Path prefix of the dataset, or the full path of a VCF file.
///
/// Returns
/// -------
/// str
///     The format's name: "VCF", "PACKEDPED" (PLINK .bed), "PED",
///     "ANCESTRYMAP", "TGENO", "PACKEDANCESTRYMAP" or "EIGENSTRAT", checked
///     in that order. No genotype files raises ``MissingFileError``.
#[pyfunction]
pub(crate) fn detect_format(prefix: &str) -> PyResult<String> {
    let format = GenoFormat::detect(prefix)?;
    Ok(format.name().to_string())
}

/// Describe a genotype dataset without reading its genotypes.
///
/// Parameters
/// ----------
/// prefix
///     Path prefix of the dataset, or the full path of a VCF file.
///
/// Returns
/// -------
/// dict
///     ``format`` ("vcf", "plink", "tgeno", "packedancestrymap" or
///     "eigenstrat"), ``n_samples`` and ``n_snps`` (int), and
///     ``populations`` (list of str, each once; sorted for a VCF, otherwise
///     in no fixed order). No genotype files raises ``MissingFileError``.
#[pyfunction]
pub(crate) fn dataset_info(prefix: &str) -> PyResult<HashMap<String, Py<PyAny>>> {
    use std::path::Path;

    Python::attach(|py| {
        let mut info: HashMap<String, Py<PyAny>> = HashMap::new();

        if let Some(path) = vcf_path(prefix) {
            let reader = open_vcf(&path)?;
            info.insert("format".to_string(), "vcf".into_py_any(py)?);
            info.insert(
                "n_samples".to_string(),
                reader.num_samples().into_py_any(py)?,
            );
            info.insert("n_snps".to_string(), reader.num_snps().into_py_any(py)?);
            // in a fixed order: sorted
            let pops: Vec<String> = reader
                .samples()
                .iter()
                .map(|s| s.population.to_string())
                .collect::<std::collections::BTreeSet<_>>()
                .into_iter()
                .collect();
            info.insert("populations".to_string(), pops.into_py_any(py)?);
            return Ok(info);
        }

        // Check for PLINK format
        if Path::new(&format!("{prefix}.bed")).exists() {
            let reader = PlinkReader::new(prefix)?;
            info.insert("format".to_string(), "plink".into_py_any(py)?);
            info.insert(
                "n_samples".to_string(),
                reader.num_samples().into_py_any(py)?,
            );
            info.insert("n_snps".to_string(), reader.num_snps().into_py_any(py)?);

            let pops: Vec<String> = reader
                .samples()
                .iter()
                .map(|s| s.population.to_string())
                .collect::<std::collections::HashSet<_>>()
                .into_iter()
                .collect();
            info.insert("populations".to_string(), pops.into_py_any(py)?);
            return Ok(info);
        }

        // Check for EIGENSTRAT/PACKEDANCESTRYMAP format
        let geno_path = format!("{prefix}.geno");
        if Path::new(&geno_path).exists() {
            if is_tgeno_format(&geno_path)? {
                let reader = TgenoReader::new(prefix)?;
                info.insert("format".to_string(), "tgeno".into_py_any(py)?);
                info.insert(
                    "n_samples".to_string(),
                    reader.samples().len().into_py_any(py)?,
                );
                info.insert("n_snps".to_string(), reader.snps().len().into_py_any(py)?);
                let pops: Vec<String> = reader
                    .samples()
                    .iter()
                    .map(|s| s.population.to_string())
                    .collect::<std::collections::HashSet<_>>()
                    .into_iter()
                    .collect();
                info.insert("populations".to_string(), pops.into_py_any(py)?);
                return Ok(info);
            }
            if is_packed_format(&geno_path)? {
                let reader = PackedAncestryMapReader::new(prefix)?;
                info.insert("format".to_string(), "packedancestrymap".into_py_any(py)?);
                info.insert(
                    "n_samples".to_string(),
                    reader.num_samples().into_py_any(py)?,
                );
                info.insert("n_snps".to_string(), reader.num_snps().into_py_any(py)?);

                let pops: Vec<String> = reader
                    .samples()
                    .iter()
                    .map(|s| s.population.to_string())
                    .collect::<std::collections::HashSet<_>>()
                    .into_iter()
                    .collect();
                info.insert("populations".to_string(), pops.into_py_any(py)?);
            } else {
                let reader = EigenstratReader::new(prefix)?;
                info.insert("format".to_string(), "eigenstrat".into_py_any(py)?);
                info.insert(
                    "n_samples".to_string(),
                    reader.num_samples().into_py_any(py)?,
                );
                info.insert("n_snps".to_string(), reader.num_snps().into_py_any(py)?);

                let pops: Vec<String> = reader
                    .samples()
                    .iter()
                    .map(|s| s.population.to_string())
                    .collect::<std::collections::HashSet<_>>()
                    .into_iter()
                    .collect();
                info.insert("populations".to_string(), pops.into_py_any(py)?);
            }
            return Ok(info);
        }

        Err(missing_file(format!(
            "No genotype files found for prefix: {prefix}"
        )))
    })
}

/// List the populations of a dataset with their numbers of samples.
///
/// Only the sample file is read (.fam, .ind, or the VCF header).
///
/// Parameters
/// ----------
/// prefix
///     Path prefix of the dataset, or the full path of a VCF file.
///
/// Returns
/// -------
/// list of tuple
///     ``(population, n_samples)`` pairs, most samples first, ties by name.
///     No genotype files raises ``MissingFileError``.
#[pyfunction]
pub(crate) fn list_populations(py: Python<'_>, prefix: &str) -> PyResult<Vec<(String, usize)>> {
    compute(py, move || {
        use std::collections::HashMap;
        use std::path::Path;

        // Helper to count populations from samples
        fn count_populations(samples: &[Sample]) -> Vec<(String, usize)> {
            let mut counts: HashMap<String, usize> = HashMap::new();
            for sample in samples {
                *counts.entry(sample.population.to_string()).or_insert(0) += 1;
            }
            let mut pop_counts: Vec<(String, usize)> = counts.into_iter().collect();
            pop_counts.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
            pop_counts
        }

        if let Some(path) = vcf_path(prefix) {
            return Ok(count_populations(open_vcf(&path)?.samples()));
        }

        // Check for PLINK format
        // only the sample files: the SNP files are not needed
        if Path::new(&format!("{prefix}.bed")).exists() {
            let samples = popmix_io::plink::read_fam_file(format!("{prefix}.fam"))?;
            return Ok(count_populations(&samples));
        }

        // Check for EIGENSTRAT/PACKEDANCESTRYMAP format
        let geno_path = format!("{prefix}.geno");
        if Path::new(&geno_path).exists() {
            if is_tgeno_format(&geno_path)? {
                let reader = TgenoReader::new(prefix)?;
                return Ok(count_populations(reader.samples()));
            }
            let samples = popmix_io::eigenstrat::read_ind_file(format!("{prefix}.ind"))?;
            return Ok(count_populations(&samples));
        }

        Err(missing_file(format!(
            "No genotype files found for prefix: {prefix}"
        )))
    })
}

/// Suggest population names of a dataset close to a mistyped one.
///
/// Matches are ranked, ignoring case: the same name in another case, then
/// names starting with the query, names containing it, names containing
/// its first three characters, and names within Levenshtein distance 3.
/// Within a rank, populations with more samples come first. The query
/// itself, when it is a population, is not suggested.
///
/// Parameters
/// ----------
/// prefix
///     Path prefix of the dataset, or the full path of a VCF file.
/// query
///     The population name to find suggestions for.
/// limit
///     Maximum number of suggestions (default 5).
///
/// Returns
/// -------
/// list of tuple
///     ``(population, n_samples)`` pairs, closest first; empty if none is
///     close.
#[pyfunction]
#[pyo3(signature = (prefix, query, *, limit=5))]
pub(crate) fn suggest_populations(
    py: Python<'_>,
    prefix: &str,
    query: &str,
    limit: usize,
) -> PyResult<Vec<(String, usize)>> {
    let all_pops = list_populations(py, prefix)?;

    if all_pops.is_empty() {
        return Ok(vec![]);
    }

    let query_lower = query.to_lowercase();

    // Score each population by similarity
    let mut scored: Vec<(i32, String, usize)> = all_pops
        .into_iter()
        .map(|(pop, count)| {
            let pop_lower = pop.to_lowercase();

            // Priority scores (lower = better match)
            let score = if pop_lower == query_lower && pop != query {
                0 // Case-different exact match (suggest the correct case)
            } else if pop_lower == query_lower {
                -1 // Exact match - skip
            } else if pop_lower.starts_with(&query_lower) {
                1 // Prefix match
            } else if pop_lower.contains(&query_lower) {
                2 // Substring match
            } else if query_lower.chars().count() >= 3
                && pop_lower.contains(&query_lower.chars().take(3).collect::<String>())
            {
                3 // First 3 chars match
            } else {
                // Levenshtein distance (higher = worse)
                let dist = levenshtein_distance(&query_lower, &pop_lower);
                if dist <= 3 {
                    4 + dist as i32 // Close typo
                } else {
                    100 + dist as i32 // Too different
                }
            };

            (score, pop, count)
        })
        .collect();

    // Sort by score, then by count descending
    scored.sort_by(|a, b| a.0.cmp(&b.0).then(b.2.cmp(&a.2)));

    // Return top matches (excluding exact matches and poor matches)
    Ok(scored
        .into_iter()
        .filter(|(score, _, _)| *score >= 0 && *score < 100)
        .take(limit)
        .map(|(_, pop, count)| (pop, count))
        .collect())
}

/// Simple Levenshtein distance implementation.
pub(crate) fn levenshtein_distance(a: &str, b: &str) -> usize {
    let a_chars: Vec<char> = a.chars().collect();
    let b_chars: Vec<char> = b.chars().collect();
    let (m, n) = (a_chars.len(), b_chars.len());

    if m == 0 {
        return n;
    }
    if n == 0 {
        return m;
    }

    let mut prev: Vec<usize> = (0..=n).collect();
    let mut curr = vec![0; n + 1];

    for i in 1..=m {
        curr[0] = i;
        for j in 1..=n {
            let cost = usize::from(a_chars[i - 1] != b_chars[j - 1]);
            curr[j] = (prev[j] + 1).min(curr[j - 1] + 1).min(prev[j - 1] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev[n]
}

// ============================================================================
// F2Blocks wrapper
// ============================================================================
