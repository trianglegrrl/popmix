//! What a dataset holds: its format, its size, its populations, and the
//! names nearest a mistyped one.

use crate::compute;
use crate::vcf_input::open_vcf;
use popmix_core::python::missing_file;
use popmix_core::Sample;
use popmix_io::tgeno::{is_tgeno_format, TgenoReader};
use popmix_io::{
    is_packed_format, EigenstratReader, GenoFormat, PackedAncestryMapReader, PlinkReader,
};
use pyo3::prelude::*;
use pyo3::IntoPyObjectExt;
use std::collections::HashMap;

/// Detect the genotype file format at a given prefix.
///
/// Args:
///     prefix: Path prefix for the dataset (or full path for VCF files)
///
/// Returns:
///     Format name: "ANCESTRYMAP", "EIGENSTRAT", "PED", "PACKEDPED", "PACKEDANCESTRYMAP", or "VCF"
#[pyfunction]
pub(crate) fn detect_format(prefix: &str) -> PyResult<String> {
    let format = GenoFormat::detect(prefix)?;
    Ok(format.name().to_string())
}

/// Get information about a genotype dataset without reading all data.
///
/// Args:
///     prefix: Path prefix for the dataset
///
/// Returns:
///     Dictionary with keys: format, n_samples, n_snps, populations
#[pyfunction]
pub(crate) fn dataset_info(prefix: &str) -> PyResult<HashMap<String, Py<PyAny>>> {
    use std::path::Path;

    Python::attach(|py| {
        let mut info: HashMap<String, Py<PyAny>> = HashMap::new();

        // Check for VCF format (full path with extension)
        if (prefix.ends_with(".vcf") || prefix.ends_with(".vcf.gz")) && Path::new(prefix).exists() {
            let reader = open_vcf(prefix)?;
            info.insert("format".to_string(), "vcf".into_py_any(py)?);
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

        // Check for VCF format with prefix
        let vcf_gz_path = format!("{}.vcf.gz", prefix);
        if Path::new(&vcf_gz_path).exists() {
            let reader = open_vcf(&vcf_gz_path)?;
            info.insert("format".to_string(), "vcf".into_py_any(py)?);
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
        let vcf_path = format!("{}.vcf", prefix);
        if Path::new(&vcf_path).exists() {
            let reader = open_vcf(&vcf_path)?;
            info.insert("format".to_string(), "vcf".into_py_any(py)?);
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

        // Check for PLINK format
        if Path::new(&format!("{}.bed", prefix)).exists() {
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
        let geno_path = format!("{}.geno", prefix);
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
            "No genotype files found for prefix: {}",
            prefix
        )))
    })
}

/// List all populations in a dataset with sample counts.
///
/// Args:
///     prefix: Path prefix for the dataset
///
/// Returns:
///     List of (population_name, sample_count) tuples, sorted by count descending
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

        // Check for VCF format (full path with extension)
        if (prefix.ends_with(".vcf") || prefix.ends_with(".vcf.gz")) && Path::new(prefix).exists() {
            let reader = open_vcf(prefix)?;
            return Ok(count_populations(reader.samples()));
        }

        // Check for VCF format with prefix
        let vcf_gz_path = format!("{}.vcf.gz", prefix);
        if Path::new(&vcf_gz_path).exists() {
            let reader = open_vcf(&vcf_gz_path)?;
            return Ok(count_populations(reader.samples()));
        }
        let vcf_path = format!("{}.vcf", prefix);
        if Path::new(&vcf_path).exists() {
            let reader = open_vcf(&vcf_path)?;
            return Ok(count_populations(reader.samples()));
        }

        // Check for PLINK format
        // only the sample files: the SNP files are not needed
        if Path::new(&format!("{}.bed", prefix)).exists() {
            let samples = popmix_io::plink::read_fam_file(&format!("{prefix}.fam"))?;
            return Ok(count_populations(&samples));
        }

        // Check for EIGENSTRAT/PACKEDANCESTRYMAP format
        let geno_path = format!("{}.geno", prefix);
        if Path::new(&geno_path).exists() {
            if is_tgeno_format(&geno_path)? {
                let reader = TgenoReader::new(prefix)?;
                return Ok(count_populations(reader.samples()));
            }
            let samples = popmix_io::eigenstrat::read_ind_file(&format!("{prefix}.ind"))?;
            return Ok(count_populations(&samples));
        }

        Err(missing_file(format!(
            "No genotype files found for prefix: {}",
            prefix
        )))
    })
}

/// Suggest similar population names for typos/errors.
///
/// Uses simple string similarity:
/// 1. Case-insensitive prefix match (highest priority)
/// 2. Case-insensitive substring match
/// 3. Levenshtein distance for typos
///
/// Args:
///     prefix: Path prefix for the dataset
///     query: The population name to find suggestions for
///     limit: Maximum number of suggestions to return (default: 5)
///
/// Returns:
///     List of (population_name, sample_count) tuples, most similar first
#[pyfunction]
#[pyo3(signature = (prefix, query, limit=5))]
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
            let cost = if a_chars[i - 1] == b_chars[j - 1] {
                0
            } else {
                1
            };
            curr[j] = (prev[j] + 1).min(curr[j - 1] + 1).min(prev[j - 1] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev[n]
}

// ============================================================================
// F2Blocks wrapper
// ============================================================================
