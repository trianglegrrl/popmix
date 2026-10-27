//! Python bindings for popmix.
//!
//! This crate provides PyO3-based Python bindings for the popmix library.

// PyO3 macro expansions generate these warnings
#![allow(clippy::useless_conversion)]

mod admixture;
mod convert;
mod geno_fstats;
mod kinship;
mod pca;
mod ploidy;
mod qpadm_models;
mod rcompat;
mod rolloff;
mod sfs;
mod vcf_input;

use std::collections::{BTreeMap, HashMap};

use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;

use convert::{convertf, merge};
use popmix_core::{AlleleFrequencies, BlockDefinition, F2Blocks, GenotypeMatrix, Ploidy, Sample};
use popmix_io::tgeno::{is_tgeno_format, TgenoReader};
use popmix_io::{
    is_packed_format, read_f2_blocks, write_f2_blocks, EigenstratReader, GenoFormat,
    PackedAncestryMapReader, PedReader, PlinkReader,
};
use vcf_input::open_vcf;

// ============================================================================
// Error conversion helper
// ============================================================================

fn to_py_err(e: popmix_core::PopmixError) -> PyErr {
    use popmix_core::PopmixError as E;
    use pyo3::exceptions::{PyFileNotFoundError, PyOSError, PyRuntimeError, PyValueError};
    match e {
        E::FileNotFound { .. } => PyFileNotFoundError::new_err(e.to_string()),
        E::Io(ref io) if io.kind() == std::io::ErrorKind::NotFound => {
            PyFileNotFoundError::new_err(e.to_string())
        }
        E::Io(_) => PyOSError::new_err(e.to_string()),
        // the caller's input: a bad name, value or file
        E::PopulationNotFound(_)
        | E::DuplicatePopulation(_)
        | E::Validation(_)
        | E::InvalidFileFormat { .. }
        | E::UnexpectedFileSize { .. }
        | E::ParseError { .. }
        | E::InsufficientPopulations { .. }
        | E::InsufficientData { .. }
        | E::DimensionMismatch { .. }
        | E::InvalidSubset { .. }
        | E::InvalidValue { .. }
        | E::InvalidGraph { .. }
        | E::InvalidAdmixtureNode { .. }
        | E::CycleDetected { .. }
        | E::InvalidConfig { .. }
        | E::InvalidThreadCount(_)
        | E::NoSnpsRemaining
        | E::NoRootNode => PyValueError::new_err(e.to_string()),
        // the computation: a singular system, no convergence, ...
        E::SingularMatrix { .. }
        | E::ConvergenceFailure { .. }
        | E::Infeasible { .. }
        | E::NumericalOverflow { .. }
        | E::Numerical(_)
        | E::LinAlg { .. } => PyRuntimeError::new_err(e.to_string()),
    }
}

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
    let have: Vec<String> = f2.populations.iter().map(|p| p.to_string()).collect();
    Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
        "not in the f2 blocks: '{}' (they have: {})",
        unknown.join("', '"),
        have.join(", ")
    )))
}

/// Requested populations each once, and each with samples in the data.
fn require_sample_populations(samples: &[popmix_core::Sample], pops: &[String]) -> PyResult<()> {
    let present: std::collections::HashSet<&str> =
        samples.iter().map(|s| s.population.as_str()).collect();
    popmix_core::genotypes::check_populations(pops.iter().map(|p| p.as_str()), |p| {
        present.contains(p)
    })
    .map_err(to_py_err)
}

// ============================================================================
// Thread pool functions
// ============================================================================

/// Set the number of threads popmix computes with (4 until set); it can be
/// changed at any time.
///
/// Args:
///     threads: Number of threads to use. Defaults to 4 if not specified.
///
/// Returns:
///     True.
#[pyfunction]
#[pyo3(signature = (threads=None))]
fn init_threads(threads: Option<usize>) -> PyResult<bool> {
    popmix_core::init_thread_pool(threads).map_err(to_py_err)
}

/// Run `f` without the GIL on popmix's thread pool (`init_threads`).
pub(crate) fn compute<T: Send>(py: Python<'_>, f: impl FnOnce() -> T + Send) -> T {
    py.detach(|| popmix_core::threading::install(f))
}

/// Get the current number of threads.
#[pyfunction]
fn get_threads() -> usize {
    popmix_core::threading::current_threads()
}

/// Get the popmix version.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

// ============================================================================
// Sample wrapper
// ============================================================================

/// A sample (individual) with metadata.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PySample {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    population: String,
    #[pyo3(get)]
    ploidy: String,
}

impl From<&popmix_core::Sample> for PySample {
    fn from(s: &popmix_core::Sample) -> Self {
        PySample {
            id: s.id.clone(),
            population: s.population.to_string(),
            ploidy: match s.ploidy {
                Ploidy::Diploid => "diploid".to_string(),
                Ploidy::Pseudohaploid => "pseudohaploid".to_string(),
            },
        }
    }
}

// ============================================================================
// SNP wrapper
// ============================================================================

/// A SNP (variant) with metadata.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PySnp {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    chromosome: String,
    #[pyo3(get)]
    position_bp: u32,
    #[pyo3(get)]
    position_cm: f64,
    #[pyo3(get)]
    allele_ref: char,
    #[pyo3(get)]
    allele_alt: char,
}

impl From<&popmix_core::Snp> for PySnp {
    fn from(s: &popmix_core::Snp) -> Self {
        PySnp {
            id: s.id.clone(),
            chromosome: s.chromosome.to_string(),
            position_bp: s.position_bp,
            // Convert from Morgans (internal storage) to cM (Python API)
            position_cm: s.position_morgans * 100.0,
            allele_ref: s.allele_ref,
            allele_alt: s.allele_alt,
        }
    }
}

// ============================================================================
// GenotypeMatrix wrapper
// ============================================================================

/// A matrix of genotypes with sample and SNP metadata.
#[pyclass]
struct PyGenotypeMatrix {
    inner: GenotypeMatrix,
}

#[pymethods]
impl PyGenotypeMatrix {
    /// Get the number of samples.
    #[getter]
    fn n_samples(&self) -> usize {
        self.inner.n_samples()
    }

    /// Get the number of SNPs.
    #[getter]
    fn n_snps(&self) -> usize {
        self.inner.n_snps()
    }

    /// Get sample metadata.
    #[getter]
    fn samples(&self) -> Vec<PySample> {
        self.inner.samples().iter().map(PySample::from).collect()
    }

    /// Get SNP metadata.
    #[getter]
    fn snps(&self) -> Vec<PySnp> {
        self.inner.snps().iter().map(PySnp::from).collect()
    }

    /// Get unique population names.
    #[getter]
    fn populations(&self) -> Vec<String> {
        self.inner
            .populations()
            .iter()
            .map(|p| p.to_string())
            .collect()
    }

    /// Get genotypes as a numpy array.
    ///
    /// Returns:
    ///     2D numpy array of shape (n_snps, n_samples) with values:
    ///     - 0 = Homozygous reference
    ///     - 1 = Heterozygous
    ///     - 2 = Homozygous alternate
    ///     - 3 = Missing
    fn to_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<u8>> {
        // one copy of the contiguous matrix
        PyArray2::from_array(py, &self.inner.data())
    }

    /// Compute allele frequencies per population.
    ///
    /// Args:
    ///     populations: Optional list of population names to compute frequencies for.
    ///                  If None, computes for all populations (slow for large datasets).
    ///
    /// Returns:
    ///     PyAlleleFrequencies with frequency data for the specified populations.
    #[pyo3(signature = (populations=None))]
    fn allele_frequencies(
        &self,
        populations: Option<Vec<String>>,
    ) -> PyResult<PyAlleleFrequencies> {
        let af = match populations {
            Some(pops) => {
                require_sample_populations(self.inner.samples(), &pops)?;
                let pop_refs: Vec<&str> = pops.iter().map(|s| s.as_str()).collect();
                AlleleFrequencies::from_genotypes_for_populations(&self.inner, Some(&pop_refs))
            }
            None => AlleleFrequencies::from_genotypes(&self.inner),
        };
        Ok(PyAlleleFrequencies { inner: af })
    }
}

// ============================================================================
// AlleleFrequencies wrapper
// ============================================================================

/// Allele frequencies computed per population.
#[pyclass]
struct PyAlleleFrequencies {
    inner: AlleleFrequencies,
}

#[pymethods]
impl PyAlleleFrequencies {
    /// Get the number of SNPs.
    #[getter]
    fn n_snps(&self) -> usize {
        self.inner.n_snps()
    }

    /// Get the number of populations.
    #[getter]
    fn n_populations(&self) -> usize {
        self.inner.n_populations()
    }

    /// Get population names.
    #[getter]
    fn populations(&self) -> Vec<String> {
        self.inner
            .populations
            .iter()
            .map(|p| p.to_string())
            .collect()
    }

    /// Get allele frequencies as a numpy array.
    ///
    /// Returns:
    ///     2D numpy array of shape (n_snps, n_populations) with frequency values.
    ///     NaN indicates missing data.
    fn to_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.inner.frequencies.clone().into_pyarray(py)
    }

    /// Get observation counts as a numpy array.
    ///
    /// Returns:
    ///     2D numpy array of shape (n_snps, n_populations) with observation counts.
    fn counts_to_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<u32>> {
        self.inner.counts.clone().into_pyarray(py)
    }

    /// Get allele frequency for a specific SNP and population.
    fn get(&self, snp_idx: usize, pop_name: &str) -> Option<f64> {
        let pop_idx = self.inner.pop_index(pop_name)?;
        self.inner.get(snp_idx, pop_idx)
    }

    /// Filter to only include specified populations.
    ///
    /// Returns a new AlleleFrequencies with only the specified populations;
    /// an unknown or repeated population is a ValueError.
    ///
    /// Args:
    ///     populations: List of population names to include
    ///
    /// Returns:
    ///     PyAlleleFrequencies with only the specified populations
    fn filter_populations(&self, populations: Vec<String>) -> PyResult<Self> {
        let pop_refs: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
        Ok(PyAlleleFrequencies {
            inner: self
                .inner
                .filter_populations(&pop_refs)
                .map_err(to_py_err)?,
        })
    }
}

// ============================================================================
// I/O functions
// ============================================================================

/// Read genotypes from PLINK binary format files.
///
/// Args:
///     prefix: Path prefix for the dataset (without .bed/.bim/.fam extension)
///     snp_range: Optional tuple (start, end) to read only a range of SNPs
///
/// Returns:
///     PyGenotypeMatrix containing the genotype data and metadata
#[pyfunction]
#[pyo3(signature = (prefix, snp_range=None))]
fn read_plink(
    py: Python<'_>,
    prefix: &str,
    snp_range: Option<(usize, usize)>,
) -> PyResult<PyGenotypeMatrix> {
    compute(py, move || {
        let reader = PlinkReader::new(prefix).map_err(to_py_err)?;

        let geno = match snp_range {
            Some((start, end)) => reader.read_snp_range(start..end).map_err(to_py_err)?,
            None => reader.read_all().map_err(to_py_err)?,
        };

        Ok(PyGenotypeMatrix { inner: geno })
    })
}

/// Read genotypes from EIGENSTRAT text format files.
///
/// Args:
///     prefix: Path prefix for the dataset (without .geno/.snp/.ind extension)
///     snp_range: Optional tuple (start, end) to read only a range of SNPs
///
/// Returns:
///     PyGenotypeMatrix containing the genotype data and metadata
#[pyfunction]
#[pyo3(signature = (prefix, snp_range=None))]
fn read_eigenstrat(
    py: Python<'_>,
    prefix: &str,
    snp_range: Option<(usize, usize)>,
) -> PyResult<PyGenotypeMatrix> {
    compute(py, move || {
        let reader = EigenstratReader::new(prefix).map_err(to_py_err)?;

        let geno = match snp_range {
            Some((start, end)) => reader.read_snp_range(start..end).map_err(to_py_err)?,
            None => reader.read_all().map_err(to_py_err)?,
        };

        Ok(PyGenotypeMatrix { inner: geno })
    })
}

/// Read genotypes from PACKEDANCESTRYMAP binary format files.
///
/// Args:
///     prefix: Path prefix for the dataset (without .geno/.snp/.ind extension)
///     snp_range: Optional tuple (start, end) to read only a range of SNPs
///
/// Returns:
///     PyGenotypeMatrix containing the genotype data and metadata
#[pyfunction]
#[pyo3(signature = (prefix, snp_range=None))]
fn read_packedancestrymap(
    py: Python<'_>,
    prefix: &str,
    snp_range: Option<(usize, usize)>,
) -> PyResult<PyGenotypeMatrix> {
    compute(py, move || {
        let reader = PackedAncestryMapReader::new(prefix).map_err(to_py_err)?;

        let geno = match snp_range {
            Some((start, end)) => reader.read_snp_range(start..end).map_err(to_py_err)?,
            None => reader.read_all().map_err(to_py_err)?,
        };

        Ok(PyGenotypeMatrix { inner: geno })
    })
}

/// Read genotypes from VCF (Variant Call Format) files.
///
/// Supports both plain text (.vcf) and gzip-compressed (.vcf.gz) files.
///
/// Args:
///     path: Path to the VCF file (including .vcf or .vcf.gz extension)
///
/// Returns:
///     PyGenotypeMatrix containing the genotype data and metadata
///
/// Note:
///     VCF files don't include population information. Each sample's
///     population is set to its sample ID by default.
#[pyfunction]
fn read_vcf(py: Python<'_>, path: &str) -> PyResult<PyGenotypeMatrix> {
    compute(py, move || {
        let reader = open_vcf(path)?;
        let geno = reader.read_all().map_err(to_py_err)?;
        Ok(PyGenotypeMatrix { inner: geno })
    })
}

/// Auto-detect file format and read genotypes.
///
/// Detects whether the dataset is PLINK, EIGENSTRAT, PACKEDANCESTRYMAP, or VCF format
/// and reads it accordingly.
///
/// Args:
///     prefix: Path prefix for the dataset
///     snp_range: Optional tuple (start, end) to read only a range of SNPs
///
/// Returns:
///     PyGenotypeMatrix containing the genotype data and metadata; samples
///     with no heterozygous call among the SNPs read are pseudohaploid
#[pyfunction]
#[pyo3(signature = (prefix, snp_range=None, populations=None))]
fn read_genotypes(
    py: Python<'_>,
    prefix: &str,
    snp_range: Option<(usize, usize)>,
    populations: Option<Vec<String>>,
) -> PyResult<PyGenotypeMatrix> {
    if let Some(pops) = populations {
        if snp_range.is_some() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "give snp_range or populations, not both",
            ));
        }
        let inner = compute(py, || read_populations(prefix, &pops))?;
        return Ok(PyGenotypeMatrix { inner });
    }
    let Some(range) = snp_range else {
        return Ok(PyGenotypeMatrix {
            inner: compute(py, || read_genotypes_internal(prefix))?,
        });
    };
    let mut geno = read_genotypes_range(py, prefix, range)?;
    ploidy::PloidyRule::Genome.apply(&mut geno.inner);
    Ok(geno)
}

/// A range of SNPs from a PLINK or (packed) EIGENSTRAT dataset; the whole
/// file for a VCF, which has no range reader.
fn read_genotypes_range(
    py: Python<'_>,
    prefix: &str,
    range: (usize, usize),
) -> PyResult<PyGenotypeMatrix> {
    use std::path::Path;
    let snp_range = Some(range);

    // Check for VCF format (full path with extension)
    if (prefix.ends_with(".vcf") || prefix.ends_with(".vcf.gz")) && Path::new(prefix).exists() {
        // VCF doesn't support snp_range, ignore it
        return read_vcf(py, prefix);
    }

    // Check for VCF format with prefix
    let vcf_gz_path = format!("{}.vcf.gz", prefix);
    if Path::new(&vcf_gz_path).exists() {
        return read_vcf(py, &vcf_gz_path);
    }
    let vcf_path = format!("{}.vcf", prefix);
    if Path::new(&vcf_path).exists() {
        return read_vcf(py, &vcf_path);
    }

    // Check for PLINK format (.bed file)
    if Path::new(&format!("{}.bed", prefix)).exists() {
        return read_plink(py, prefix, snp_range);
    }

    // Check for EIGENSTRAT/PACKEDANCESTRYMAP format (.geno file)
    let geno_path = format!("{}.geno", prefix);
    if Path::new(&geno_path).exists() {
        // Check if it's packed (binary), transposed packed or text format
        if is_tgeno_format(&geno_path).map_err(to_py_err)? {
            let reader = TgenoReader::new(prefix).map_err(to_py_err)?;
            let (start, end) = range;
            return Ok(PyGenotypeMatrix {
                inner: reader.read_snp_range(start..end).map_err(to_py_err)?,
            });
        }
        if is_packed_format(&geno_path).map_err(to_py_err)? {
            return read_packedancestrymap(py, prefix, snp_range);
        } else {
            return read_eigenstrat(py, prefix, snp_range);
        }
    }

    Err(PyErr::new::<pyo3::exceptions::PyFileNotFoundError, _>(
        format!("No genotype files found for prefix: {}", prefix),
    ))
}

/// Load precomputed F2 blocks from a directory.
///
/// The directory should contain files created by F2Blocks.save():
/// - metadata.json: Population names and block info
/// - block_sums.npy: NumPy array of F2 sums
/// - block_weights.npy: NumPy array of weights
///
/// Args:
///     path: Directory path containing the F2 block files
///
/// Returns:
///     F2Blocks object ready for computing f-statistics
#[pyfunction]
fn read_f2(py: Python<'_>, path: &str) -> PyResult<PyF2Blocks> {
    compute(py, move || {
        let inner = read_f2_blocks(std::path::Path::new(path)).map_err(to_py_err)?;
        Ok(PyF2Blocks { inner })
    })
}

/// Detect the genotype file format at a given prefix.
///
/// Args:
///     prefix: Path prefix for the dataset (or full path for VCF files)
///
/// Returns:
///     Format name: "ANCESTRYMAP", "EIGENSTRAT", "PED", "PACKEDPED", "PACKEDANCESTRYMAP", or "VCF"
#[pyfunction]
fn detect_format(prefix: &str) -> PyResult<String> {
    let format = GenoFormat::detect(prefix).map_err(to_py_err)?;
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
fn dataset_info(prefix: &str) -> PyResult<HashMap<String, Py<PyAny>>> {
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
            let reader = PlinkReader::new(prefix).map_err(to_py_err)?;
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
            if is_tgeno_format(&geno_path).map_err(to_py_err)? {
                let reader = TgenoReader::new(prefix).map_err(to_py_err)?;
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
            if is_packed_format(&geno_path).map_err(to_py_err)? {
                let reader = PackedAncestryMapReader::new(prefix).map_err(to_py_err)?;
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
                let reader = EigenstratReader::new(prefix).map_err(to_py_err)?;
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

        Err(PyErr::new::<pyo3::exceptions::PyFileNotFoundError, _>(
            format!("No genotype files found for prefix: {}", prefix),
        ))
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
fn list_populations(py: Python<'_>, prefix: &str) -> PyResult<Vec<(String, usize)>> {
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
            let samples =
                popmix_io::plink::read_fam_file(&format!("{prefix}.fam")).map_err(to_py_err)?;
            return Ok(count_populations(&samples));
        }

        // Check for EIGENSTRAT/PACKEDANCESTRYMAP format
        let geno_path = format!("{}.geno", prefix);
        if Path::new(&geno_path).exists() {
            if is_tgeno_format(&geno_path).map_err(to_py_err)? {
                let reader = TgenoReader::new(prefix).map_err(to_py_err)?;
                return Ok(count_populations(reader.samples()));
            }
            let samples = popmix_io::eigenstrat::read_ind_file(&format!("{prefix}.ind"))
                .map_err(to_py_err)?;
            return Ok(count_populations(&samples));
        }

        Err(PyErr::new::<pyo3::exceptions::PyFileNotFoundError, _>(
            format!("No genotype files found for prefix: {}", prefix),
        ))
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
fn suggest_populations(
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
fn levenshtein_distance(a: &str, b: &str) -> usize {
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

/// Precomputed F2 statistics for all population pairs.
///
/// This is the main data structure for computing f-statistics.
/// Create using `extract_f2()` or from allele frequencies.
#[pyclass]
struct PyF2Blocks {
    inner: F2Blocks,
}

#[pymethods]
impl PyF2Blocks {
    /// Get the number of populations.
    #[getter]
    fn n_populations(&self) -> usize {
        self.inner.n_populations()
    }

    /// Get the number of blocks.
    #[getter]
    fn n_blocks(&self) -> usize {
        self.inner.n_blocks
    }

    /// These populations only, in this order (admixtools2
    /// ``extract_f2_subset``); FST blocks are kept.
    fn subset(&self, populations: Vec<String>) -> PyResult<PyF2Blocks> {
        let refs: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
        Ok(PyF2Blocks {
            inner: self.inner.subset_populations(&refs).map_err(to_py_err)?,
        })
    }

    /// SNPs with data for both populations in each block (the counts
    /// admixtools2 stores with its f2 blocks), or None for unknown
    /// populations.
    fn counts_per_block(&self, pop_a: &str, pop_b: &str) -> Option<Vec<f64>> {
        self.inner.block_counts(pop_a, pop_b)
    }

    /// F2 blocks from each pair's per-block f2 mean and SNP count, as
    /// admixtools2 stores them: ``means[i][j]`` and ``counts[i][j]`` are the
    /// blocks of populations i and j (only i <= j is read). FST is not
    /// available from these.
    #[staticmethod]
    fn from_block_means(
        populations: Vec<String>,
        block_lengths: Vec<usize>,
        means: Vec<Vec<Vec<f64>>>,
        counts: Vec<Vec<Vec<f64>>>,
    ) -> PyResult<PyF2Blocks> {
        let pops = populations
            .iter()
            .map(|p| popmix_core::PopulationId::from(p.as_str()))
            .collect();
        Ok(PyF2Blocks {
            inner: F2Blocks::from_block_means(pops, block_lengths, &means, &counts)
                .map_err(to_py_err)?,
        })
    }

    /// These blocks, in this order and with repeats, as admixtools2
    /// ``f2_blocks[, , blocks]`` (0-based). FST blocks are not carried over.
    fn select_blocks(&self, blocks: Vec<usize>) -> PyResult<PyF2Blocks> {
        Ok(PyF2Blocks {
            inner: self.inner.select_blocks(&blocks).map_err(to_py_err)?,
        })
    }

    /// Number of SNPs in each block (the jackknife weights).
    #[getter]
    fn block_lengths(&self) -> Vec<usize> {
        self.inner.block_lengths.clone()
    }

    /// Number of SNPs in each FST block (FST uses its own SNP set and blocks,
    /// as admixtools2). Empty if the blocks carry no FST data (f2 directories
    /// saved before storage version 4).
    #[getter]
    fn fst_block_lengths(&self) -> Vec<usize> {
        self.inner
            .fst
            .as_ref()
            .map(|f| f.block_lengths.clone())
            .unwrap_or_default()
    }

    /// Per-block FST for a pair of populations, or None if unavailable.
    fn fst_per_block<'py>(
        &self,
        py: Python<'py>,
        pop_a: &str,
        pop_b: &str,
    ) -> Option<Bound<'py, PyArray1<f64>>> {
        let fst = self.inner.fst.as_ref()?;
        let pair = self
            .inner
            .pair_index(self.inner.pop_index(pop_a)?, self.inner.pop_index(pop_b)?);
        Some(ndarray::Array1::from(fst.block_values(pair)).into_pyarray(py))
    }

    /// Get population names.
    #[getter]
    fn populations(&self) -> Vec<String> {
        self.inner
            .populations
            .iter()
            .map(|p| p.to_string())
            .collect()
    }

    /// Compute f2(A, B) statistic.
    ///
    /// Args:
    ///     pop_a: First population name
    ///     pop_b: Second population name
    ///
    /// Returns:
    ///     Dictionary with keys: estimate, std_err, z_score, n_snps
    fn f2(&self, pop_a: &str, pop_b: &str) -> PyResult<Option<BTreeMap<String, f64>>> {
        require_populations(&self.inner, &[pop_a, pop_b])?;
        match self.inner.f2(pop_a, pop_b) {
            Some((estimate, std_err, n_snps)) => {
                let mut result = BTreeMap::new();
                result.insert("estimate".to_string(), estimate);
                result.insert("std_err".to_string(), std_err);
                result.insert(
                    "z_score".to_string(),
                    if std_err > 0.0 {
                        estimate / std_err
                    } else {
                        0.0
                    },
                );
                result.insert("n_snps".to_string(), n_snps as f64);
                Ok(Some(result))
            }
            None => Ok(None),
        }
    }

    /// Compute f3(A; B, C) statistic (admixture test).
    ///
    /// A negative f3 indicates that A is admixed between B and C.
    ///
    /// Args:
    ///     pop_a: Target population (potentially admixed)
    ///     pop_b: Source population 1
    ///     pop_c: Source population 2
    ///
    /// Returns:
    ///     Dictionary with keys: estimate, std_err, z_score, n_blocks
    fn f3(&self, pop_a: &str, pop_b: &str, pop_c: &str) -> PyResult<Option<BTreeMap<String, f64>>> {
        require_populations(&self.inner, &[pop_a, pop_b, pop_c])?;
        match self.inner.f3(pop_a, pop_b, pop_c) {
            Some(result) => {
                let mut out = BTreeMap::new();
                out.insert("estimate".to_string(), result.estimate);
                out.insert("std_err".to_string(), result.std_err);
                out.insert("z_score".to_string(), result.z_score);
                out.insert("n_blocks".to_string(), result.n_blocks as f64);
                Ok(Some(out))
            }
            None => Ok(None),
        }
    }

    /// Compute f4(A, B; C, D) statistic (D-statistic).
    ///
    /// Tests for gene flow between (A,B) and (C,D) branches.
    ///
    /// Args:
    ///     pop_a, pop_b: First pair of populations
    ///     pop_c, pop_d: Second pair of populations
    ///
    /// Returns:
    ///     Dictionary with keys: estimate, std_err, z_score, n_blocks
    fn f4(
        &self,
        pop_a: &str,
        pop_b: &str,
        pop_c: &str,
        pop_d: &str,
    ) -> PyResult<Option<BTreeMap<String, f64>>> {
        require_populations(&self.inner, &[pop_a, pop_b, pop_c, pop_d])?;
        match self.inner.f4(pop_a, pop_b, pop_c, pop_d) {
            Some(result) => {
                let mut out = BTreeMap::new();
                out.insert("estimate".to_string(), result.estimate);
                out.insert("std_err".to_string(), result.std_err);
                out.insert("z_score".to_string(), result.z_score);
                out.insert("n_blocks".to_string(), result.n_blocks as f64);
                Ok(Some(out))
            }
            None => Ok(None),
        }
    }

    /// f4(A, B; C, D) from the f2 blocks, under its old name.
    ///
    /// This is f4, not Patterson's D: D is normalized per SNP and needs the
    /// genotypes (``popmix.qpdstat(prefix, ..., f4mode=False)``).
    fn d_stat(
        &self,
        pop_a: &str,
        pop_b: &str,
        pop_c: &str,
        pop_d: &str,
    ) -> PyResult<Option<BTreeMap<String, f64>>> {
        self.f4(pop_a, pop_b, pop_c, pop_d)
    }

    /// Compute FST (Hudson's fixation index) between two populations.
    ///
    /// FST measures population differentiation on a scale of 0 to 1:
    /// - FST = 0: No differentiation (populations are identical)
    /// - FST = 1: Complete fixation (populations share no alleles)
    ///
    /// Args:
    ///     pop_a: First population name
    ///     pop_b: Second population name
    ///
    /// Returns:
    ///     Dictionary with keys: estimate, std_err, n_snps
    fn fst(&self, pop_a: &str, pop_b: &str) -> PyResult<Option<BTreeMap<String, f64>>> {
        require_populations(&self.inner, &[pop_a, pop_b])?;
        match self.inner.fst(pop_a, pop_b) {
            Some(result) => {
                let mut out = BTreeMap::new();
                out.insert("estimate".to_string(), result.estimate);
                out.insert("std_err".to_string(), result.std_err);
                out.insert("n_snps".to_string(), result.n_snps as f64);
                Ok(Some(out))
            }
            None => Ok(None),
        }
    }

    /// Compute D-statistics (qpdstat) for multiple population combinations.
    ///
    /// This tests for gene flow between population pairs. D(A, B; C, D) is
    /// positive if there's more gene flow from D to A (or C to B).
    ///
    /// Args:
    ///     pop1: First population set (list of population names)
    ///     pop2: Second population set (list of population names)
    ///     pop3: Third population set (list of population names)
    ///     pop4: Fourth population set (list of population names)
    ///
    /// Returns:
    ///     List of dictionaries with keys: pop1, pop2, pop3, pop4, estimate, std_err, z_score
    fn qpdstat(
        &self,
        pop1: Vec<String>,
        pop2: Vec<String>,
        pop3: Vec<String>,
        pop4: Vec<String>,
    ) -> PyResult<Vec<HashMap<String, Py<PyAny>>>> {
        let pop1_refs: Vec<&str> = pop1.iter().map(|s| s.as_str()).collect();
        let pop2_refs: Vec<&str> = pop2.iter().map(|s| s.as_str()).collect();
        let pop3_refs: Vec<&str> = pop3.iter().map(|s| s.as_str()).collect();
        let pop4_refs: Vec<&str> = pop4.iter().map(|s| s.as_str()).collect();
        Python::attach(|py| {
            self.inner
                .qpdstat(&pop1_refs, &pop2_refs, &pop3_refs, &pop4_refs)
                .into_iter()
                .map(|(p1, p2, p3, p4, result)| {
                    let mut out = HashMap::new();
                    out.insert("pop1".to_string(), p1.into_py_any(py)?);
                    out.insert("pop2".to_string(), p2.into_py_any(py)?);
                    out.insert("pop3".to_string(), p3.into_py_any(py)?);
                    out.insert("pop4".to_string(), p4.into_py_any(py)?);
                    out.insert("estimate".to_string(), result.estimate.into_py_any(py)?);
                    out.insert("std_err".to_string(), result.std_err.into_py_any(py)?);
                    out.insert("z_score".to_string(), result.z_score.into_py_any(py)?);
                    out.insert(
                        "n_blocks".to_string(),
                        (result.n_blocks as f64).into_py_any(py)?,
                    );
                    Ok(out)
                })
                .collect()
        })
    }

    /// Compute f3 outgroup statistics (qp3pop) for multiple population combinations.
    ///
    /// f3(Target; Source1, Source2) tests if Target is admixed between Source1
    /// and Source2. A significantly negative f3 indicates admixture.
    ///
    /// Args:
    ///     targets: Target populations to test for admixture
    ///     source1: First source population set
    ///     source2: Second source population set
    ///
    /// Returns:
    ///     List of dictionaries with keys: target, source1, source2, estimate, std_err, z_score
    fn qp3pop(
        &self,
        targets: Vec<String>,
        source1: Vec<String>,
        source2: Vec<String>,
    ) -> PyResult<Vec<HashMap<String, Py<PyAny>>>> {
        let targets_refs: Vec<&str> = targets.iter().map(|s| s.as_str()).collect();
        let source1_refs: Vec<&str> = source1.iter().map(|s| s.as_str()).collect();
        let source2_refs: Vec<&str> = source2.iter().map(|s| s.as_str()).collect();
        Python::attach(|py| {
            self.inner
                .qp3pop(&targets_refs, &source1_refs, &source2_refs)
                .into_iter()
                .map(|(target, s1, s2, result)| {
                    let mut out = HashMap::new();
                    out.insert("target".to_string(), target.into_py_any(py)?);
                    out.insert("source1".to_string(), s1.into_py_any(py)?);
                    out.insert("source2".to_string(), s2.into_py_any(py)?);
                    out.insert("estimate".to_string(), result.estimate.into_py_any(py)?);
                    out.insert("std_err".to_string(), result.std_err.into_py_any(py)?);
                    out.insert("z_score".to_string(), result.z_score.into_py_any(py)?);
                    out.insert(
                        "n_blocks".to_string(),
                        (result.n_blocks as f64).into_py_any(py)?,
                    );
                    Ok(out)
                })
                .collect()
        })
    }

    /// Get the f2 matrix as a numpy array.
    ///
    /// Returns:
    ///     2D numpy array of shape (n_pops, n_pops) with f2 values.
    fn f2_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.inner.f2_matrix().into_pyarray(py)
    }

    /// Get per-block f2 values for a pair of populations.
    ///
    /// Returns:
    ///     1D numpy array of shape (n_blocks,) with f2 values for each block.
    fn f2_per_block<'py>(
        &self,
        py: Python<'py>,
        pop_a: &str,
        pop_b: &str,
    ) -> PyResult<Option<Bound<'py, PyArray1<f64>>>> {
        match self.inner.get_block_f2_values(pop_a, pop_b) {
            Some(values) => Ok(Some(values.into_pyarray(py))),
            None => Ok(None),
        }
    }

    /// Save F2 blocks to a directory.
    ///
    /// Creates a directory with:
    /// - metadata.json: Population names, block info
    /// - block_sums.npy: NumPy format array of F2 sums
    /// - block_weights.npy: NumPy format array of weights
    ///
    /// Args:
    ///     path: Directory path to write to (will be created if needed)
    ///     source: Optional source file path for metadata
    #[pyo3(signature = (path, source=None))]
    fn save(&self, path: &str, source: Option<&str>) -> PyResult<()> {
        write_f2_blocks(&self.inner, std::path::Path::new(path), source).map_err(to_py_err)
    }
}

// ============================================================================
// F-statistics functions
// ============================================================================

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
fn extract_f2(
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
    compute(py, move || {
        use std::path::Path;

        let relabel = sample_labels(inds, pops)?;
        let rule = ploidy::ploidy_rule(pseudohaploid)?;

        if !(0.0..=1.0).contains(&maxmiss) {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "maxmiss must be between 0 and 1, got {maxmiss}"
            )));
        }
        let f2_options = popmix_core::F2Options {
            max_missing: maxmiss,
            poly_only,
        };

        // Create block definition
        let block_def = match block_size {
            Some(size) => BlockDefinition::new(size),
            None => BlockDefinition::default(),
        };

        // Packed .geno and PLINK .bed are read in place, SNP by SNP
        let geno_path = format!("{}.geno", prefix);
        if let Some(ploidy::Mapped {
            file,
            mut samples,
            snps,
        }) = ploidy::mapped(prefix)?
        {
            popmix_core::genotypes::check_sorted(&snps).map_err(to_py_err)?;
            if let Some(labels) = &relabel {
                apply_sample_labels(&mut samples, labels)?;
            }
            require_sample_populations(&samples, &populations)?;
            let pop_refs: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
            rule.apply_mapped(&file, &mut samples, &pop_refs);
            // a chromosome at a time: only one chromosome's frequency table
            // is held (blocks never span chromosomes)
            let mut runs: Vec<std::ops::Range<usize>> = Vec::new();
            for (i, snp) in snps.iter().enumerate() {
                match runs.last_mut() {
                    Some(r) if snps[r.start].chromosome == snp.chromosome => r.end = i + 1,
                    _ => runs.push(i..i + 1),
                }
            }
            let mut evidence = popmix_core::CoordinateEvidence::default();
            for snp in snps.iter().filter(|s| !auto_only || s.is_autosomal()) {
                evidence.observe(std::slice::from_ref(snp));
            }
            let parts = runs.iter().map(|r| {
                file.allele_frequencies_at(
                    &samples,
                    &snps[r.clone()],
                    r.start,
                    &pop_refs,
                    auto_only,
                )
            });
            let f2_blocks = match F2Blocks::from_parts(parts, &evidence, &block_def, &f2_options) {
                Some(blocks) => blocks,
                None => {
                    let af = file.allele_frequencies(&samples, &snps, &pop_refs, auto_only);
                    F2Blocks::from_allele_frequencies_with_options(&af, &block_def, &f2_options)
                }
            };
            return Ok(PyF2Blocks { inner: f2_blocks });
        }

        let mut geno =
            if Path::new(&geno_path).exists() && is_tgeno_format(&geno_path).map_err(to_py_err)? {
                // TGENO stores each sample's genotypes together: read only the
                // records of the requested populations' samples
                let reader = TgenoReader::new(prefix).map_err(to_py_err)?;
                let mut samples = reader.samples().to_vec();
                if let Some(labels) = &relabel {
                    apply_sample_labels(&mut samples, labels)?;
                }
                require_sample_populations(&samples, &populations)?;
                let selected: Vec<usize> = (0..samples.len())
                    .filter(|&i| {
                        populations
                            .iter()
                            .any(|p| p == samples[i].population.as_str())
                    })
                    .collect();
                let mut geno = reader.read_samples(&selected).map_err(to_py_err)?;
                for (sample, &i) in geno.samples_mut().iter_mut().zip(&selected) {
                    *sample = samples[i].clone();
                }
                geno
            } else {
                // Fallback: load all genotypes for other formats (PLINK, text EIGENSTRAT)
                let mut geno = read_genotypes_internal(prefix)?;
                if let Some(labels) = &relabel {
                    apply_sample_labels(geno.samples_mut(), labels)?;
                }
                require_sample_populations(geno.samples(), &populations)?;
                geno
            };
        popmix_core::genotypes::check_sorted(geno.snps()).map_err(to_py_err)?;
        rule.apply(&mut geno);

        // Filter to autosomal if requested
        let geno = if auto_only {
            geno.filter_autosomal().map_err(to_py_err)?
        } else {
            geno
        };

        // Compute allele frequencies for specified populations
        let pop_refs: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
        let af = AlleleFrequencies::from_genotypes_for_populations(&geno, Some(&pop_refs));

        // Compute F2 blocks
        let f2_blocks =
            F2Blocks::from_allele_frequencies_with_options(&af, &block_def, &f2_options);

        Ok(PyF2Blocks { inner: f2_blocks })
    })
}

/// The population label of an excluded sample: no population is called this.
const EXCLUDED: &str = "\u{0}excluded";

/// admixtools2 `inds`/`pops`: each listed sample's population.
fn sample_labels(
    inds: Option<Vec<String>>,
    pops: Option<Vec<String>>,
) -> PyResult<Option<HashMap<String, String>>> {
    match (inds, pops) {
        (None, None) => Ok(None),
        (Some(inds), Some(pops)) if inds.len() == pops.len() => {
            let n = inds.len();
            let labels: HashMap<String, String> = inds.into_iter().zip(pops).collect();
            if labels.len() < n {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "a sample is listed twice in inds",
                ));
            }
            Ok(Some(labels))
        }
        (Some(_), Some(_)) => Err(pyo3::exceptions::PyValueError::new_err(
            "inds and pops must have the same length: one population per sample",
        )),
        _ => Err(pyo3::exceptions::PyValueError::new_err(
            "give inds and pops together",
        )),
    }
}

/// Give the listed samples their populations and leave out the rest.
fn apply_sample_labels(samples: &mut [Sample], labels: &HashMap<String, String>) -> PyResult<()> {
    let known: std::collections::HashSet<&str> = samples.iter().map(|s| s.id.as_str()).collect();
    if let Some(missing) = labels.keys().find(|id| !known.contains(id.as_str())) {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "sample {missing:?} is not in the genotype data"
        )));
    }
    for s in samples.iter_mut() {
        s.population =
            popmix_core::PopulationId::from(labels.get(&s.id).map_or(EXCLUDED, |p| p.as_str()));
    }
    Ok(())
}

/// Compute F2 blocks from allele frequencies.
#[pyfunction]
#[pyo3(signature = (af, block_size=None))]
fn f2_from_frequencies(
    py: Python<'_>,
    af: &PyAlleleFrequencies,
    block_size: Option<f64>,
) -> PyResult<PyF2Blocks> {
    compute(py, move || {
        let block_def = match block_size {
            Some(size) => BlockDefinition::new(size),
            None => BlockDefinition::default(),
        };

        let f2_blocks = F2Blocks::from_allele_frequencies(&af.inner, &block_def);

        Ok(PyF2Blocks { inner: f2_blocks })
    })
}

// ============================================================================
// qpWave/qpAdm wrappers
// ============================================================================

/// Result of a single rank test from qpWave.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyRankTest {
    #[pyo3(get)]
    rank: usize,
    #[pyo3(get)]
    chi_sq: f64,
    #[pyo3(get)]
    df: usize,
    #[pyo3(get)]
    p_value: f64,
}

/// Result of qpWave analysis.
#[pyclass]
struct PyQpWaveResult {
    #[pyo3(get)]
    rank_tests: Vec<PyRankTest>,
    #[pyo3(get)]
    estimated_rank: usize,
    #[pyo3(get)]
    left_pops: Vec<String>,
    #[pyo3(get)]
    right_pops: Vec<String>,
    /// f4(left[i+1], left[0]; right[j+1], right[0]), as rows i, columns j.
    #[pyo3(get)]
    f4_matrix: Vec<Vec<f64>>,
    /// Jackknife covariance of the f4 matrix flattened by rows (left first).
    #[pyo3(get)]
    f4_covariance: Vec<Vec<f64>>,
}

#[pymethods]
impl PyQpWaveResult {
    /// Check if populations form a clade (rank 0 p-value > threshold).
    #[pyo3(signature = (p_threshold=0.05))]
    fn is_clade(&self, p_threshold: f64) -> bool {
        self.rank_tests
            .first()
            .map(|r| r.p_value > p_threshold)
            .unwrap_or(false)
    }

    /// Get the full rank test.
    fn full_rank(&self) -> Option<PyRankTest> {
        self.rank_tests.last().cloned()
    }
}

/// Result of f4 ratio analysis for admixture proportion estimation.
#[pyclass]
struct PyF4RatioResult {
    #[pyo3(get)]
    test: String,
    #[pyo3(get)]
    source_a: String,
    #[pyo3(get)]
    source_b: String,
    #[pyo3(get)]
    outgroup: String,
    #[pyo3(get)]
    reference: String,
    #[pyo3(get)]
    alpha: f64,
    #[pyo3(get)]
    std_err: f64,
    #[pyo3(get)]
    z_score: f64,
    #[pyo3(get)]
    n_blocks: usize,
    #[pyo3(get)]
    f4_numerator: f64,
    #[pyo3(get)]
    f4_denominator: f64,
}

#[pymethods]
impl PyF4RatioResult {
    fn __repr__(&self) -> String {
        format!(
            "F4RatioResult(test='{}', source_a='{}', source_b='{}', alpha={:.4}, std_err={:.4})",
            self.test, self.source_a, self.source_b, self.alpha, self.std_err
        )
    }
}

/// Estimate admixture proportions via f4 ratios.
///
/// Computes alpha = f4(A,O;X,C) / f4(A,O;B,C) where:
/// - X is the test (admixed) population
/// - A is source population contributing alpha proportion
/// - B is source population contributing (1-alpha) proportion
/// - O is an outgroup
/// - C is a reference population
///
/// Args:
///     f2_blocks: Pre-computed F2 blocks
///     test: Test (admixed) population X
///     source_a: Source population A
///     source_b: Source population B
///     outgroup: Outgroup population O
///     reference: Reference population C
///
/// Returns:
///     PyF4RatioResult with estimated alpha and standard error
#[pyfunction]
#[pyo3(signature = (f2_blocks, test, source_a, source_b, outgroup, reference))]
fn qpf4ratio(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    test: String,
    source_a: String,
    source_b: String,
    outgroup: String,
    reference: String,
) -> PyResult<PyF4RatioResult> {
    compute(py, move || {
        let result = popmix_core::f4ratio(
            &f2_blocks.inner,
            &test,
            &source_a,
            &source_b,
            &outgroup,
            &reference,
            None,
        )
        .map_err(to_py_err)?;

        Ok(PyF4RatioResult {
            test: result.test.to_string(),
            source_a: result.source_a.to_string(),
            source_b: result.source_b.to_string(),
            outgroup: result.outgroup.to_string(),
            reference: result.reference.to_string(),
            alpha: result.alpha,
            std_err: result.std_err,
            z_score: result.z_score,
            n_blocks: result.n_blocks,
            f4_numerator: result.f4_numerator,
            f4_denominator: result.f4_denominator,
        })
    })
}

/// Result of qpAdm analysis.
#[pyclass]
struct PyQpAdmResult {
    #[pyo3(get)]
    target: String,
    #[pyo3(get)]
    sources: Vec<String>,
    #[pyo3(get)]
    weights: Vec<f64>,
    #[pyo3(get)]
    weight_std_errors: Vec<f64>,
    #[pyo3(get)]
    chi_sq: f64,
    #[pyo3(get)]
    df: usize,
    #[pyo3(get)]
    p_value: f64,
    #[pyo3(get)]
    constrained: bool,
    /// Fit of every subset of the sources (admixtools2 ``popdrop``).
    #[pyo3(get)]
    popdrop: Vec<PyPopDropRow>,
}

/// One row of qpAdm ``popdrop``: the model fitted to a subset of the sources.
///
/// ``pattern`` has one character per source, "1" if it is dropped; ``weights``
/// is None for dropped sources. ``best`` marks the nested path from the full
/// model, each step dropping one more source; ``dofdiff``, ``chisqdiff`` and
/// ``p_nested`` compare a step with the model it came from. With
/// ``compat="admixtools2-2.0.10"`` they reproduce admixtools2's columns
/// instead (DEVIATIONS.md, ADMIXTOOLS2-0001).
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyPopDropRow {
    #[pyo3(get)]
    pattern: String,
    #[pyo3(get)]
    wt: usize,
    #[pyo3(get)]
    dof: usize,
    #[pyo3(get)]
    chi_sq: f64,
    #[pyo3(get)]
    p_value: f64,
    #[pyo3(get)]
    f4rank: usize,
    #[pyo3(get)]
    weights: Vec<Option<f64>>,
    #[pyo3(get)]
    feasible: bool,
    #[pyo3(get)]
    best: bool,
    #[pyo3(get)]
    dofdiff: Option<usize>,
    #[pyo3(get)]
    chisqdiff: Option<f64>,
    #[pyo3(get)]
    p_nested: Option<f64>,
}

impl From<popmix_core::QpAdmResult> for PyQpAdmResult {
    fn from(result: popmix_core::QpAdmResult) -> Self {
        PyQpAdmResult {
            target: result.target.to_string(),
            sources: result.sources.iter().map(|s| s.to_string()).collect(),
            weights: result.weights,
            weight_std_errors: result.weight_std_errors,
            chi_sq: result.chi_sq,
            df: result.df,
            p_value: result.p_value,
            constrained: result.constrained,
            popdrop: result
                .popdrop
                .into_iter()
                .map(|r| PyPopDropRow {
                    pattern: r.pattern,
                    wt: r.wt,
                    dof: r.dof,
                    chi_sq: r.chi_sq,
                    p_value: r.p_value,
                    f4rank: r.f4rank,
                    weights: r.weights,
                    feasible: r.feasible,
                    best: r.best,
                    dofdiff: r.dofdiff,
                    chisqdiff: r.chisqdiff,
                    p_nested: r.p_nested,
                })
                .collect(),
        }
    }
}

#[pymethods]
impl PyQpAdmResult {
    /// Check if the model is a good fit (p > threshold).
    #[pyo3(signature = (p_threshold=0.05))]
    fn is_good_fit(&self, p_threshold: f64) -> bool {
        self.p_value > p_threshold
    }

    /// Check if all weights are feasible (in [0, 1]).
    fn is_feasible(&self) -> bool {
        self.weights.iter().all(|&w| (0.0..=1.0).contains(&w))
    }

    /// Get weights as a dictionary mapping source name to weight, in source order.
    fn weights_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        for (source, weight) in self.sources.iter().zip(&self.weights) {
            dict.set_item(source, weight)?;
        }
        Ok(dict)
    }
}

/// A single model result from qpAdm rotation.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyRotationModel {
    #[pyo3(get)]
    sources: Vec<String>,
    #[pyo3(get)]
    weights: Vec<f64>,
    #[pyo3(get)]
    weight_std_errors: Vec<f64>,
    #[pyo3(get)]
    p_value: f64,
    #[pyo3(get)]
    chi_sq: f64,
    #[pyo3(get)]
    df: usize,
    #[pyo3(get)]
    is_feasible: bool,
}

#[pymethods]
impl PyRotationModel {
    fn __repr__(&self) -> String {
        format!(
            "RotationModel(sources={:?}, p_value={:.4}, is_feasible={})",
            self.sources, self.p_value, self.is_feasible
        )
    }

    /// Get weights as a dictionary mapping source name to weight, in source order.
    fn weights_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        for (source, weight) in self.sources.iter().zip(&self.weights) {
            dict.set_item(source, weight)?;
        }
        Ok(dict)
    }
}

/// Result of qpAdm rotation analysis.
#[pyclass]
struct PyQpAdmRotateResult {
    #[pyo3(get)]
    target: String,
    #[pyo3(get)]
    candidates: Vec<String>,
    #[pyo3(get)]
    outgroups: Vec<String>,
    #[pyo3(get)]
    models: Vec<PyRotationModel>,
    #[pyo3(get)]
    n_combinations: usize,
    #[pyo3(get)]
    n_feasible: usize,
    /// Combinations not fit: more sources than outgroups, or a singular or
    /// infeasible system.
    #[pyo3(get)]
    n_skipped: usize,
}

#[pymethods]
impl PyQpAdmRotateResult {
    /// Get the best model (highest p-value).
    fn best_model(&self) -> Option<PyRotationModel> {
        self.models.first().cloned()
    }

    /// Get models passing a p-value threshold.
    #[pyo3(signature = (p_threshold=0.05))]
    fn passing_models(&self, p_threshold: f64) -> Vec<PyRotationModel> {
        self.models
            .iter()
            .filter(|m| m.p_value > p_threshold)
            .cloned()
            .collect()
    }

    /// Get feasible models only.
    fn feasible_models(&self) -> Vec<PyRotationModel> {
        self.models
            .iter()
            .filter(|m| m.is_feasible)
            .cloned()
            .collect()
    }
}

/// Run qpAdm rotation to test all source combinations.
///
/// Tests all combinations of candidate source populations to find
/// the best-fitting admixture models.
///
/// Args:
///     f2_blocks: Precomputed F2 blocks
///     target: Target population name
///     candidates: Candidate source population names
///     outgroups: Outgroup population names
///     min_sources: Minimum number of sources to test (default: 1)
///     max_sources: Maximum number of sources to test (default: all)
///     feasible_only: Only return feasible models (default: False)
///     fudge: Regularization factor (default: 1e-4)
///
/// Returns:
///     PyQpAdmRotateResult with all tested models sorted by p-value
#[pyfunction]
#[pyo3(signature = (f2_blocks, target, candidates, outgroups, min_sources=None, max_sources=None, feasible_only=None, fudge=None))]
#[allow(clippy::too_many_arguments)]
fn qpadm_rotate(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    target: String,
    candidates: Vec<String>,
    outgroups: Vec<String>,
    min_sources: Option<usize>,
    max_sources: Option<usize>,
    feasible_only: Option<bool>,
    fudge: Option<f64>,
) -> PyResult<PyQpAdmRotateResult> {
    compute(py, move || {
        let candidates_refs: Vec<&str> = candidates.iter().map(|s| s.as_str()).collect();
        let outgroups_refs: Vec<&str> = outgroups.iter().map(|s| s.as_str()).collect();

        let opts = popmix_core::QpAdmRotateOptions {
            min_sources: min_sources.unwrap_or(1),
            max_sources,
            p_threshold: 0.05,
            feasible_only: feasible_only.unwrap_or(false),
            qpadm_opts: popmix_core::QpAdmOptions {
                fudge: fudge.unwrap_or(1e-4),
                ..Default::default()
            },
        };

        let result = popmix_core::qpadm_rotate(
            &f2_blocks.inner,
            &target,
            &candidates_refs,
            &outgroups_refs,
            Some(opts),
        )
        .map_err(to_py_err)?;

        let models = result
            .models
            .into_iter()
            .map(|m| PyRotationModel {
                sources: m.sources.iter().map(|s| s.to_string()).collect(),
                weights: m.weights,
                weight_std_errors: m.weight_std_errors,
                p_value: m.p_value,
                chi_sq: m.chi_sq,
                df: m.df,
                is_feasible: m.is_feasible,
            })
            .collect();

        Ok(PyQpAdmRotateResult {
            target: result.target.to_string(),
            candidates: result.candidates.iter().map(|s| s.to_string()).collect(),
            outgroups: result.outgroups.iter().map(|s| s.to_string()).collect(),
            models,
            n_combinations: result.n_combinations,
            n_feasible: result.n_feasible,
            n_skipped: result.n_skipped,
        })
    })
}

/// Run qpWave analysis.
///
/// Tests whether left populations can be modeled as mixtures of ancestral
/// populations distinguished by the right (outgroup) populations.
///
/// Args:
///     f2_blocks: Precomputed F2 blocks
///     left: Left population names
///     right: Right (outgroup) population names
///     fudge: Regularization factor (default: 1e-4)
///
/// Returns:
///     PyQpWaveResult with rank tests and estimated rank
#[pyfunction]
#[pyo3(signature = (f2_blocks, left, right, fudge=None))]
fn qpwave(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    left: Vec<String>,
    right: Vec<String>,
    fudge: Option<f64>,
) -> PyResult<PyQpWaveResult> {
    compute(py, move || {
        let left_refs: Vec<&str> = left.iter().map(|s| s.as_str()).collect();
        let right_refs: Vec<&str> = right.iter().map(|s| s.as_str()).collect();

        let opts = popmix_core::QpWaveOptions {
            fudge: fudge.unwrap_or(1e-4),
            max_rank: None,
        };

        let result = popmix_core::qpwave(&f2_blocks.inner, &left_refs, &right_refs, Some(opts))
            .map_err(to_py_err)?;

        let rank_tests = result
            .rank_tests
            .iter()
            .map(|r| PyRankTest {
                rank: r.rank,
                chi_sq: r.chi_sq,
                df: r.df,
                p_value: r.p_value,
            })
            .collect();

        let rows = |m: &ndarray::Array2<f64>| m.outer_iter().map(|r| r.to_vec()).collect();
        Ok(PyQpWaveResult {
            rank_tests,
            estimated_rank: result.estimated_rank,
            left_pops: left,
            right_pops: right,
            f4_matrix: rows(&result.f4_matrix.values),
            f4_covariance: rows(&result.f4_matrix.covariance),
        })
    })
}

/// Result of qpWave test for a single population pair.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyQpWavePairResult {
    #[pyo3(get)]
    pop_a: String,
    #[pyo3(get)]
    pop_b: String,
    #[pyo3(get)]
    p_value: f64,
    #[pyo3(get)]
    chi_sq: f64,
    #[pyo3(get)]
    df: usize,
}

#[pymethods]
impl PyQpWavePairResult {
    fn __repr__(&self) -> String {
        format!(
            "QpWavePairResult(pop_a='{}', pop_b='{}', p_value={:.4})",
            self.pop_a, self.pop_b, self.p_value
        )
    }

    /// Check if this pair forms a clade (p > threshold).
    #[pyo3(signature = (p_threshold=0.05))]
    fn is_clade(&self, p_threshold: f64) -> bool {
        self.p_value > p_threshold
    }
}

/// Result of pairwise qpWave analysis.
#[pyclass]
struct PyQpWavePairsResult {
    #[pyo3(get)]
    populations: Vec<String>,
    #[pyo3(get)]
    outgroups: Vec<String>,
    #[pyo3(get)]
    pairs: Vec<PyQpWavePairResult>,
    #[pyo3(get)]
    n_pairs: usize,
    /// Pairs whose qpWave system was singular, left out of `pairs`.
    #[pyo3(get)]
    n_skipped: usize,
}

#[pymethods]
impl PyQpWavePairsResult {
    /// Get pairs that pass the cladality test (p > threshold).
    #[pyo3(signature = (p_threshold=0.05))]
    fn clades(&self, p_threshold: f64) -> Vec<PyQpWavePairResult> {
        self.pairs
            .iter()
            .filter(|p| p.p_value > p_threshold)
            .cloned()
            .collect()
    }

    /// Get the pair with highest p-value (most likely clade).
    fn best_clade(&self) -> Option<PyQpWavePairResult> {
        self.pairs
            .iter()
            .filter(|p| !p.p_value.is_nan())
            .max_by(|a, b| a.p_value.total_cmp(&b.p_value))
            .cloned()
    }
}

/// Run qpWave for all population pairs (cladality screening).
///
/// Tests whether each pair of populations forms a clade relative to
/// the outgroup populations.
///
/// Args:
///     f2_blocks: Precomputed F2 blocks
///     populations: Population names to test pairwise
///     outgroups: Outgroup population names
///     fudge: Regularization factor (default: 1e-4)
///
/// Returns:
///     PyQpWavePairsResult with all pairwise test results
#[pyfunction]
#[pyo3(signature = (f2_blocks, populations, outgroups, fudge=None))]
fn qpwave_pairs(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    populations: Vec<String>,
    outgroups: Vec<String>,
    fudge: Option<f64>,
) -> PyResult<PyQpWavePairsResult> {
    compute(py, move || {
        let pop_refs: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
        let outgroup_refs: Vec<&str> = outgroups.iter().map(|s| s.as_str()).collect();

        let opts = popmix_core::QpWaveOptions {
            fudge: fudge.unwrap_or(1e-4),
            max_rank: None,
        };

        let result =
            popmix_core::qpwave_pairs(&f2_blocks.inner, &pop_refs, &outgroup_refs, Some(opts))
                .map_err(to_py_err)?;

        let pairs = result
            .pairs
            .into_iter()
            .map(|p| PyQpWavePairResult {
                pop_a: p.pop_a.to_string(),
                pop_b: p.pop_b.to_string(),
                p_value: p.p_value,
                chi_sq: p.chi_sq,
                df: p.df,
            })
            .collect();

        Ok(PyQpWavePairsResult {
            populations: result.populations.iter().map(|s| s.to_string()).collect(),
            outgroups: result.outgroups.iter().map(|s| s.to_string()).collect(),
            pairs,
            n_pairs: result.n_pairs,
            n_skipped: result.n_skipped,
        })
    })
}

/// Parse a qpAdm `compat` argument.
fn admixtools2_compat(compat: Option<&str>) -> PyResult<popmix_core::Admixtools2Compat> {
    match compat {
        None => Ok(popmix_core::Admixtools2Compat::CORRECTED),
        Some("admixtools2-2.0.10") => Ok(popmix_core::Admixtools2Compat::V2_0_10),
        Some(other) => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unknown compat {other:?} for qpadm; expected None or \"admixtools2-2.0.10\""
        ))),
    }
}

/// Run qpAdm analysis to estimate admixture weights.
///
/// Estimates the proportion of ancestry in a target population from
/// source populations, using outgroup populations to constrain the model.
///
/// Args:
///     f2_blocks: Precomputed F2 blocks
///     target: Target population name
///     sources: Source population names
///     right: Right (outgroup) population names
///     constrained: Whether to constrain weights to [0, 1] (default: False)
///     fudge: Regularization factor (default: 1e-4)
///     compat: Reproduce admixtools2's known quirks instead of correcting
///         them: "admixtools2-2.0.10", or None (default) for popmix's
///         corrected behaviour. See DEVIATIONS.md.
///
/// Returns:
///     PyQpAdmResult with estimated weights and model fit statistics
#[pyfunction]
#[pyo3(signature = (f2_blocks, target, sources, right, constrained=None, fudge=None, compat=None))]
#[allow(clippy::too_many_arguments)]
fn qpadm(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    target: String,
    sources: Vec<String>,
    right: Vec<String>,
    constrained: Option<bool>,
    fudge: Option<f64>,
    compat: Option<&str>,
) -> PyResult<PyQpAdmResult> {
    compute(py, move || {
        let sources_refs: Vec<&str> = sources.iter().map(|s| s.as_str()).collect();
        let right_refs: Vec<&str> = right.iter().map(|s| s.as_str()).collect();

        let opts = popmix_core::QpAdmOptions {
            fudge: fudge.unwrap_or(1e-4),
            constrained: constrained.unwrap_or(false),
            compat: admixtools2_compat(compat)?,
        };

        let result = popmix_core::qpadm(
            &f2_blocks.inner,
            &target,
            &sources_refs,
            &right_refs,
            Some(opts),
        )
        .map_err(to_py_err)?;

        Ok(PyQpAdmResult::from(result))
    })
}

/// Run qpAdm for multiple targets with the same sources and outgroups.
///
/// This is a convenience function for testing multiple target populations
/// against the same admixture model.
///
/// Args:
///     f2_blocks: Pre-computed F2 blocks
///     targets: List of target populations to test
///     sources: Source populations for the admixture model
///     right: Right (outgroup) populations
///     constrained: Use constrained optimization (weights in [0,1])
///     fudge: Regularization parameter
///     compat: "admixtools2-2.0.10" to reproduce admixtools2's quirks (see
///         qpadm)
///
/// Returns:
///     List of qpAdm results for each target
#[pyfunction]
#[pyo3(signature = (f2_blocks, targets, sources, right, constrained=None, fudge=None, compat=None))]
#[allow(clippy::too_many_arguments)]
fn qpadm_multi(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    targets: Vec<String>,
    sources: Vec<String>,
    right: Vec<String>,
    constrained: Option<bool>,
    fudge: Option<f64>,
    compat: Option<&str>,
) -> PyResult<Vec<PyQpAdmResult>> {
    compute(py, move || {
        let sources_refs: Vec<&str> = sources.iter().map(|s| s.as_str()).collect();
        let right_refs: Vec<&str> = right.iter().map(|s| s.as_str()).collect();

        let opts = popmix_core::QpAdmOptions {
            fudge: fudge.unwrap_or(1e-4),
            constrained: constrained.unwrap_or(false),
            compat: admixtools2_compat(compat)?,
        };

        let mut results = Vec::with_capacity(targets.len());

        for target in &targets {
            match popmix_core::qpadm(
                &f2_blocks.inner,
                target,
                &sources_refs,
                &right_refs,
                Some(opts.clone()),
            ) {
                Ok(result) => {
                    results.push(PyQpAdmResult::from(result));
                }
                Err(e) => {
                    return Err(to_py_err(e));
                }
            }
        }

        Ok(results)
    })
}

// ============================================================================
// qpGraph wrappers
// ============================================================================

/// An admixture graph representing population relationships.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyAdmixGraph {
    inner: popmix_core::AdmixGraph,
}

#[pymethods]
impl PyAdmixGraph {
    /// Get the number of nodes.
    #[getter]
    fn node_count(&self) -> usize {
        self.inner.node_count()
    }

    /// Get the number of edges.
    #[getter]
    fn edge_count(&self) -> usize {
        self.inner.edge_count()
    }

    /// Get leaf population names.
    #[getter]
    fn leaves(&self) -> Vec<String> {
        self.inner.leaf_names()
    }

    /// Get number of drift edges.
    #[getter]
    fn num_drift_edges(&self) -> usize {
        self.inner.num_drift_edges()
    }

    /// Get edge weights.
    fn get_edge_weights(&self) -> Vec<f64> {
        self.inner.get_edge_weights()
    }

    /// Compute expected f2 between two populations.
    fn expected_f2(&self, pop_a: &str, pop_b: &str) -> PyResult<f64> {
        self.inner.expected_f2(pop_a, pop_b).map_err(to_py_err)
    }
}

/// Build a simple tree graph using the builder pattern.
///
/// Args:
///     root_name: Name of the root node
///     leaves: List of (name, parent, drift) tuples for leaf nodes
///     internals: Optional list of (name, parent, drift) tuples for internal nodes
///
/// Returns:
///     PyAdmixGraph
#[pyfunction]
#[pyo3(signature = (root_name, leaves, internals=None))]
fn build_tree(
    py: Python<'_>,
    root_name: &str,
    leaves: Vec<(String, String, f64)>,
    internals: Option<Vec<(String, String, f64)>>,
) -> PyResult<PyAdmixGraph> {
    compute(py, move || {
        let mut builder = popmix_core::GraphBuilder::new().root(root_name);

        // Add internal nodes first
        if let Some(ints) = internals {
            for (name, parent, drift) in ints {
                builder = builder.internal(&name, &parent, drift);
            }
        }

        // Add leaves
        for (name, parent, drift) in leaves {
            builder = builder.leaf(&name, &parent, drift);
        }

        let graph = builder.build_unchecked();

        Ok(PyAdmixGraph { inner: graph })
    })
}

/// One graph found by find_graphs.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyScoredGraph {
    /// The graph as (from, to) pairs; the path from the root to the first
    /// population comes first, so it is the f3 base.
    #[pyo3(get)]
    edges: Vec<(String, String)>,
    /// qpgraph score.
    #[pyo3(get)]
    score: f64,
}

#[pymethods]
impl PyScoredGraph {
    fn __repr__(&self) -> String {
        format!(
            "ScoredGraph(score={:.6}, {} edges)",
            self.score,
            self.edges.len()
        )
    }
}

/// Search for the admixture graphs that best fit the data.
///
/// Every candidate graph has the populations as leaves and ``nadmix``
/// admixture nodes and is fitted with qpgraph(); ``populations[0]`` is the f3
/// base population of every fit.
///
/// Args:
///     f2_blocks: Precomputed F2 blocks containing every population
///     populations: The graph's leaves
///     nadmix: Number of admixture nodes (default: 0)
///     outgroup: Only graphs whose root has this population as a direct child
///     method: "stepwise" (default): add the populations, then the admixture
///         edges, one at a time, keeping the ``beam_width`` best graphs after
///         each step. "exhaustive": every graph admixtools2
///         generate_all_graphs() returns (feasible for about 5 populations
///         with 1 admixture, or 4 with 2).
///     beam_width: Graphs kept after each stepwise step (default: 20)
///     numstart, seed, diag, diag_f3, constrained: as for qpgraph()
///
/// Returns:
///     (graphs, n_evaluated): the graphs of the last step, best first, and
///     the number of qpgraph fits made
#[pyfunction]
#[pyo3(signature = (f2_blocks, populations, nadmix=0, outgroup=None, method="stepwise", beam_width=20, numstart=10, seed=1, diag=1e-4, diag_f3=1e-5, constrained=true))]
#[allow(clippy::too_many_arguments)]
fn find_graphs(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    populations: Vec<String>,
    nadmix: usize,
    outgroup: Option<String>,
    method: &str,
    beam_width: usize,
    numstart: usize,
    seed: i32,
    diag: f64,
    diag_f3: f64,
    constrained: bool,
) -> PyResult<(Vec<PyScoredGraph>, usize)> {
    compute(py, move || {
        let method = match method {
            "stepwise" => popmix_core::GraphSearchMethod::Stepwise { beam_width },
            "exhaustive" => popmix_core::GraphSearchMethod::Exhaustive,
            other => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "method must be \"stepwise\" or \"exhaustive\", got {other:?}"
                )))
            }
        };
        let opts = popmix_core::GraphSearchOptions {
            nadmix,
            outgroup,
            method,
            qpgraph: popmix_core::QpGraphOptions {
                diag,
                diag_f3,
                constrained,
                numstart,
                seed,
                admix_weights: None,
                f2_blocks_test: None,
            },
        };
        let pop_refs: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
        let result =
            popmix_core::find_graphs(&f2_blocks.inner, &pop_refs, Some(opts)).map_err(to_py_err)?;
        let graphs = result
            .graphs
            .into_iter()
            .map(|g| PyScoredGraph {
                edges: g.edges,
                score: g.score,
            })
            .collect();
        Ok((graphs, result.n_evaluated))
    })
}

/// Every admixture graph with these leaves and ``nadmix`` admixture nodes,
/// as admixtools2 generate_all_graphs() returns them, in its order.
///
/// Each graph is a list of (from, to) pairs whose path from the root to
/// ``leaves[0]`` comes first.
#[pyfunction]
#[pyo3(signature = (leaves, nadmix=0))]
fn all_graphs(
    py: Python<'_>,
    leaves: Vec<String>,
    nadmix: usize,
) -> PyResult<Vec<Vec<(String, String)>>> {
    compute(py, move || {
        let refs: Vec<&str> = leaves.iter().map(|s| s.as_str()).collect();
        popmix_core::all_graphs(&refs, nadmix).map_err(to_py_err)
    })
}

/// The value admixtools2 graph_hash() digests for a simplified graph: two
/// graphs are the same graph to admixtools2 exactly when their keys are equal.
#[pyfunction]
fn graph_key(edges: Vec<(String, String)>) -> PyResult<Vec<String>> {
    popmix_core::graph_key(&edges).map_err(to_py_err)
}

/// One fitted qpgraph edge.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyQpGraphEdge {
    #[pyo3(get)]
    from_node: String,
    #[pyo3(get)]
    to_node: String,
    /// "edge" (drift) or "admix", as in admixtools2.
    #[pyo3(get)]
    kind: String,
    /// Drift length, or admixture weight.
    #[pyo3(get)]
    weight: f64,
    /// Admixture edges: lowest and highest weight over the random starts.
    #[pyo3(get)]
    low: Option<f64>,
    #[pyo3(get)]
    high: Option<f64>,
}

/// One fitted f3(pop1; pop2, pop3) of a qpgraph fit.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyQpGraphF3 {
    #[pyo3(get)]
    pop1: String,
    #[pyo3(get)]
    pop2: String,
    #[pyo3(get)]
    pop3: String,
    #[pyo3(get)]
    est: f64,
    #[pyo3(get)]
    se: f64,
    #[pyo3(get)]
    fit: f64,
    #[pyo3(get)]
    diff: f64,
    #[pyo3(get)]
    z: f64,
    #[pyo3(get)]
    p: f64,
}

/// One random start of the qpgraph admixture weight optimization.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyQpGraphStart {
    #[pyo3(get)]
    start: Vec<f64>,
    #[pyo3(get)]
    weights: Vec<f64>,
    #[pyo3(get)]
    score: f64,
    /// Score evaluations (admixtools2 ``fevals``)
    #[pyo3(get)]
    evaluations: usize,
    /// R ``optim`` convergence code: 0 converged, 1 iteration limit
    #[pyo3(get)]
    convergence: i32,
}

#[pymethods]
impl PyQpGraphStart {
    /// Whether this start converged (``convergence == 0``)
    #[getter]
    fn converged(&self) -> bool {
        self.convergence == 0
    }
}

/// Result of qpgraph().
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyQpGraphResult {
    /// Leaves in graph order; the first is the f3 base population.
    #[pyo3(get)]
    pops: Vec<String>,
    #[pyo3(get)]
    score: f64,
    #[pyo3(get)]
    edges: Vec<PyQpGraphEdge>,
    #[pyo3(get)]
    f3: Vec<PyQpGraphF3>,
    #[pyo3(get)]
    opt: Vec<PyQpGraphStart>,
    /// Inverse of the regularized f3 covariance, rows in ``f3`` order.
    #[pyo3(get)]
    ppinv: Vec<Vec<f64>>,
    /// Score against ``f2_blocks_test``, when given.
    #[pyo3(get)]
    score_test: Option<f64>,
}

#[pymethods]
impl PyQpGraphResult {
    fn __repr__(&self) -> String {
        format!(
            "QpGraphResult(score={:.4}, pops={}, edges={})",
            self.score,
            self.pops.len(),
            self.edges.len()
        )
    }
}

/// Fit an admixture graph, as admixtools2 ``qpgraph`` does.
///
/// Drift lengths minimise the covariance-weighted residual of f3(base; i, j)
/// for all pairs of leaves, with lengths >= 0; admixture weights are optimized
/// over [0, 1] from ``numstart`` random starts. The populations are the
/// graph's leaves in order of first appearance in ``edges``; the first is the
/// f3 base population.
///
/// Args:
///     f2_blocks: Precomputed F2 blocks containing every leaf
///     edges: The graph as (from, to) pairs. A node with two incoming edges is
///         an admixture node.
///     admix_weights: Fit at these weights instead of optimizing them, one
///         per admixture node in order of first appearance: the weight of the
///         incoming edge whose source node appears first.
///     numstart: Random starts for the admixture weights (default: 10)
///     seed: Seed for the random starts, as R's ``set.seed`` (default: 1);
///         with the same seed the starts are admixtools2's
///     diag: Regularization of the edge-length fit (default: 1e-4)
///     diag_f3: Regularization of the f3 covariance (default: 1e-5)
///     constrained: Lengths >= 0 and weights in [0, 1] (default: True)
///     f2_blocks_test: Also score the fit against these blocks' f3
///         statistics (``score_test``), for an out-of-sample score
///
/// Returns:
///     PyQpGraphResult with fitted edges, score, f3 fits and the random starts
#[pyfunction]
#[pyo3(signature = (f2_blocks, edges, admix_weights=None, numstart=10, seed=1, diag=1e-4, diag_f3=1e-5, constrained=true, f2_blocks_test=None))]
#[allow(clippy::too_many_arguments)]
fn qpgraph(
    py: Python<'_>,
    f2_blocks: &PyF2Blocks,
    edges: Vec<(String, String)>,
    admix_weights: Option<Vec<f64>>,
    numstart: usize,
    seed: i32,
    diag: f64,
    diag_f3: f64,
    constrained: bool,
    f2_blocks_test: Option<&PyF2Blocks>,
) -> PyResult<PyQpGraphResult> {
    compute(py, move || {
        let opts = popmix_core::QpGraphOptions {
            diag,
            diag_f3,
            constrained,
            numstart,
            seed,
            admix_weights,
            f2_blocks_test: f2_blocks_test.map(|b| b.inner.clone()),
        };
        let r = popmix_core::qpgraph(&f2_blocks.inner, &edges, Some(opts)).map_err(to_py_err)?;
        Ok(PyQpGraphResult {
            pops: r.pops,
            score: r.score,
            edges: r
                .edges
                .into_iter()
                .map(|e| PyQpGraphEdge {
                    from_node: e.from,
                    to_node: e.to,
                    kind: if e.admix { "admix" } else { "edge" }.to_string(),
                    weight: e.weight,
                    low: e.low,
                    high: e.high,
                })
                .collect(),
            f3: r
                .f3
                .into_iter()
                .map(|f| PyQpGraphF3 {
                    pop1: f.pop1,
                    pop2: f.pop2,
                    pop3: f.pop3,
                    est: f.est,
                    se: f.se,
                    fit: f.fit,
                    diff: f.diff,
                    z: f.z,
                    p: f.p,
                })
                .collect(),
            opt: r
                .opt
                .into_iter()
                .map(|o| PyQpGraphStart {
                    start: o.start,
                    weights: o.weights,
                    score: o.score,
                    evaluations: o.evaluations,
                    convergence: o.convergence,
                })
                .collect(),
            ppinv: r.ppinv.outer_iter().map(|row| row.to_vec()).collect(),
            score_test: r.score_test,
        })
    })
}

/// Every genotype of a dataset in any format, each sample pseudohaploid
/// when it has no heterozygous call (the default `pseudohaploid` rule).
/// Every SNP of the samples of `populations`, in file order: packed .geno
/// and PLINK .bed are read in place, other formats whole and then subset.
fn read_populations(prefix: &str, populations: &[String]) -> PyResult<GenotypeMatrix> {
    let (mut source, snps, _, samples) = geno_fstats::open_source(prefix, populations, false)?;
    let mut data = ndarray::Array2::<u8>::zeros((snps.len(), samples.len()));
    let mut rows = geno_fstats::rows_of(&mut source);
    for mut row in data.rows_mut() {
        rows.next_row(row.as_slice_mut().expect("row-major"))
            .map_err(to_py_err)?;
    }
    drop(rows);
    GenotypeMatrix::new(data, samples, snps).map_err(to_py_err)
}

fn read_genotypes_internal(prefix: &str) -> PyResult<GenotypeMatrix> {
    let mut geno = read_genotypes_any_format(prefix)?;
    ploidy::PloidyRule::Genome.apply(&mut geno);
    Ok(geno)
}

/// The genotypes of the samples of `populations` (or more): for TGENO only
/// their records are read, for other formats the whole file.
pub(crate) fn read_genotypes_of(prefix: &str, populations: &[String]) -> PyResult<GenotypeMatrix> {
    let geno_path = format!("{prefix}.geno");
    if std::path::Path::new(&geno_path).exists()
        && is_tgeno_format(&geno_path).map_err(to_py_err)?
    {
        let reader = TgenoReader::new(prefix).map_err(to_py_err)?;
        let selected: Vec<usize> = (0..reader.samples().len())
            .filter(|&i| {
                let pop = reader.samples()[i].population.as_str();
                populations.iter().any(|p| p == pop)
            })
            .collect();
        return reader.read_samples(&selected).map_err(to_py_err);
    }
    read_genotypes_any_format(prefix)
}

pub(crate) fn read_genotypes_any_format(prefix: &str) -> PyResult<GenotypeMatrix> {
    use std::path::Path;

    // Check for VCF format (full path with extension)
    if (prefix.ends_with(".vcf") || prefix.ends_with(".vcf.gz")) && Path::new(prefix).exists() {
        let reader = open_vcf(prefix)?;
        return reader.read_all().map_err(to_py_err);
    }

    // Check for VCF format with prefix
    let vcf_gz_path = format!("{}.vcf.gz", prefix);
    if Path::new(&vcf_gz_path).exists() {
        let reader = open_vcf(&vcf_gz_path)?;
        return reader.read_all().map_err(to_py_err);
    }
    let vcf_path = format!("{}.vcf", prefix);
    if Path::new(&vcf_path).exists() {
        let reader = open_vcf(&vcf_path)?;
        return reader.read_all().map_err(to_py_err);
    }

    // Check for PLINK format (.bed file)
    if Path::new(&format!("{}.bed", prefix)).exists() {
        let reader = PlinkReader::new(prefix).map_err(to_py_err)?;
        return reader.read_all().map_err(to_py_err);
    }

    // Text PLINK (.ped and .map)
    if Path::new(&format!("{}.ped", prefix)).exists() {
        let reader = PedReader::new(prefix).map_err(to_py_err)?;
        return reader.read_all().map_err(to_py_err);
    }

    // Check for EIGENSTRAT/PACKEDANCESTRYMAP format (.geno file)
    let geno_path = format!("{}.geno", prefix);
    if Path::new(&geno_path).exists() {
        if is_tgeno_format(&geno_path).map_err(to_py_err)? {
            let reader = TgenoReader::new(prefix).map_err(to_py_err)?;
            return reader.read_all().map_err(to_py_err);
        }
        if is_packed_format(&geno_path).map_err(to_py_err)? {
            let reader = PackedAncestryMapReader::new(prefix).map_err(to_py_err)?;
            return reader.read_all().map_err(to_py_err);
        } else {
            let reader = EigenstratReader::new(prefix).map_err(to_py_err)?;
            return reader.read_all().map_err(to_py_err);
        }
    }

    Err(PyErr::new::<pyo3::exceptions::PyFileNotFoundError, _>(
        format!("No genotype files found for prefix: {}", prefix),
    ))
}

// ============================================================================
// DATES wrappers
// ============================================================================

/// Result of a single LD bin from DATES.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyLdBin {
    #[pyo3(get)]
    distance_cm: f64,
    #[pyo3(get)]
    weighted_ld: f64,
    #[pyo3(get)]
    z_score: f64,
    #[pyo3(get)]
    correlation: f64,
    #[pyo3(get)]
    regression_coef: f64,
    #[pyo3(get)]
    n_pairs: usize,
    #[pyo3(get)]
    std_err: f64,
}

impl From<&popmix_core::LdBin> for PyLdBin {
    fn from(b: &popmix_core::LdBin) -> Self {
        PyLdBin {
            distance_cm: b.distance_cm,
            weighted_ld: b.weighted_ld,
            z_score: b.z_score,
            correlation: b.correlation,
            regression_coef: b.regression_coef,
            n_pairs: b.n_pairs,
            std_err: b.std_err,
        }
    }
}

/// Result of exponential decay fitting.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
struct PyDecayFit {
    #[pyo3(get)]
    amplitude: f64,
    #[pyo3(get)]
    decay_rate: f64,
    #[pyo3(get)]
    intercept: f64,
    #[pyo3(get)]
    r_squared: f64,
    #[pyo3(get)]
    decay_rate_se: f64,
    #[pyo3(get)]
    rss: f64,
}

impl From<&popmix_core::DecayFit> for PyDecayFit {
    fn from(f: &popmix_core::DecayFit) -> Self {
        PyDecayFit {
            amplitude: f.amplitude,
            decay_rate: f.decay_rate,
            intercept: f.intercept,
            r_squared: f.r_squared,
            decay_rate_se: f.decay_rate_se,
            rss: f.rss,
        }
    }
}

/// Result of DATES admixture dating analysis.
#[pyclass]
struct PyDatesResult {
    #[pyo3(get)]
    test_pop: String,
    #[pyo3(get)]
    source1: String,
    #[pyo3(get)]
    source2: String,
    #[pyo3(get)]
    generations: f64,
    #[pyo3(get)]
    generations_se: f64,
    #[pyo3(get)]
    years_bp: f64,
    #[pyo3(get)]
    years_bp_se: f64,
    #[pyo3(get)]
    n_snps: usize,
    #[pyo3(get)]
    n_chromosomes: usize,
    ld_bins: Vec<popmix_core::LdBin>,
    fit: popmix_core::DecayFit,
    jackknife: Vec<popmix_core::dates::JackknifeBlock>,
}

#[pymethods]
impl PyDatesResult {
    /// Per-chromosome jackknife blocks as (chromosome, weight, generations),
    /// matching the DATES `.jin` file. Empty if the jackknife did not run.
    #[getter]
    fn jackknife_by_chromosome(&self) -> Vec<(u8, f64, f64)> {
        self.jackknife
            .iter()
            .map(|b| (b.chromosome, b.weight, b.generations))
            .collect()
    }

    /// Get the LD decay bins.
    #[getter]
    fn ld_bins(&self) -> Vec<PyLdBin> {
        self.ld_bins.iter().map(PyLdBin::from).collect()
    }

    /// Get the exponential fit result.
    #[getter]
    fn fit(&self) -> PyDecayFit {
        PyDecayFit::from(&self.fit)
    }

    /// Get LD decay data as a list of dictionaries for easy pandas conversion.
    fn to_ld_dataframe<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        self.ld_bins
            .iter()
            .map(|b| {
                let d = PyDict::new(py);
                d.set_item("distance_cm", b.distance_cm)?;
                d.set_item("weighted_ld", b.weighted_ld)?;
                d.set_item("z_score", b.z_score)?;
                d.set_item("correlation", b.correlation)?;
                d.set_item("regression_coef", b.regression_coef)?;
                d.set_item("n_pairs", b.n_pairs as f64)?;
                d.set_item("std_err", b.std_err)?;
                Ok(d)
            })
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "DatesResult(test_pop='{}', generations={:.1}±{:.1}, years_bp={:.0}±{:.0}, n_snps={})",
            self.test_pop,
            self.generations,
            self.generations_se,
            self.years_bp,
            self.years_bp_se,
            self.n_snps
        )
    }
}

/// Estimate admixture time using the DATES algorithm.
///
/// DATES estimates when admixture occurred by measuring the decay of
/// linkage disequilibrium with genetic distance. LD from recent
/// admixture decays exponentially, with the decay rate proportional
/// to the number of generations since admixture.
///
/// Args:
///     geno: Genotype matrix with all relevant populations
///     test_pop: Name of the test (admixed) population
///     source1: Name of first source population
///     source2: Name of second source population
///     generation_time: Generation time in years (default: 29)
///     min_distance_cm: Minimum distance for fitting, DATES lovalfit (default: 0.45)
///     max_distance_cm: Maximum distance, DATES maxdis x 100 (default: 100.0)
///     bin_width_cm: Width of distance bins, DATES binsize x 100 (default: 0.1)
///     run_jackknife: Run jackknife for SE estimation (default: True)
///     qbin: DATES qbin. 0 visits every SNP pair; > 0 uses an FFT mesh of
///         qbin cells per bin, which is faster and slightly approximate
///         (default: 0, the DATES default)
///     compat: Reproduce a reference implementation's known quirks instead
///         of correcting them. None (default) corrects them all;
///         "dates-v753" gives output identical to DATES v753, for checking
///         against it. See DEVIATIONS.md for what each quirk changes.
///
/// Returns:
///     DatesResult with estimated admixture timing. ``generations`` and
///     ``generations_se`` are the weighted jackknife estimate and SE that
///     DATES reports; ``fit.decay_rate`` is the full-data fit.
///
/// Example:
///     >>> geno = pm.read_packedancestrymap("/path/to/data")
///     >>> result = pm.dates(geno, "French", "Yamnaya", "Anatolia_N")
///     >>> print(f"Admixture time: {result.generations:.0f} generations")
///     >>> print(f"Years BP: {result.years_bp:.0f}")
#[pyfunction]
#[pyo3(signature = (geno, test_pop, source1, source2, generation_time=None, min_distance_cm=None, max_distance_cm=None, bin_width_cm=None, run_jackknife=None, qbin=None, compat=None))]
#[allow(clippy::too_many_arguments)]
fn dates(
    py: Python<'_>,
    geno: &PyGenotypeMatrix,
    test_pop: &str,
    source1: &str,
    source2: &str,
    generation_time: Option<f64>,
    min_distance_cm: Option<f64>,
    max_distance_cm: Option<f64>,
    bin_width_cm: Option<f64>,
    run_jackknife: Option<bool>,
    qbin: Option<usize>,
    compat: Option<&str>,
) -> PyResult<PyDatesResult> {
    compute(py, move || {
        let compat = match compat {
            None => popmix_core::dates::DatesCompat::CORRECTED,
            Some("dates-v753") => popmix_core::dates::DatesCompat::V753,
            Some(other) => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "unknown compat {other:?} for dates; expected None or \"dates-v753\""
                )))
            }
        };
        let mut ld_options = popmix_core::LdOptions::default();
        if let Some(q) = qbin {
            ld_options.qbin = q;
        }
        if let Some(bw) = bin_width_cm {
            ld_options.bin_width_cm = bw;
        }
        if let Some(max_d) = max_distance_cm {
            ld_options.max_distance_cm = max_d;
        }

        let mut fit_options = popmix_core::dates::FitOptions::default();
        if let Some(bw) = bin_width_cm {
            fit_options.bin_size_cm = bw;
        }
        if let Some(min_d) = min_distance_cm {
            fit_options.min_distance_cm = min_d;
        }
        if let Some(max_d) = max_distance_cm {
            fit_options.max_distance_cm = max_d;
        }

        let opts = popmix_core::DatesOptions {
            ld_options,
            fit_options,
            generation_time: generation_time.unwrap_or(29.0),
            run_jackknife: run_jackknife.unwrap_or(true),
            compat,
            ..Default::default()
        };

        let result = popmix_core::dates(&geno.inner, None, test_pop, source1, source2, Some(opts))
            .map_err(to_py_err)?;

        Ok(PyDatesResult {
            test_pop: result.test_pop.to_string(),
            source1: result.source1.to_string(),
            source2: result.source2.to_string(),
            generations: result.generations,
            generations_se: result.generations_se,
            years_bp: result.years_bp,
            years_bp_se: result.years_bp_se,
            n_snps: result.n_snps,
            n_chromosomes: result.n_chromosomes,
            ld_bins: result.ld_bins,
            fit: result.fit,
            jackknife: result.jackknife_by_chromosome,
        })
    })
}

// ============================================================================
// Module definition
// ============================================================================

/// Python module for popmix population genetics toolkit.
// Not yet tested on free-threaded Python (3.13t): until it is, importing
// popmix there turns the GIL back on rather than risk a race.
#[pymodule(gil_used = true)]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Threading functions
    m.add_function(wrap_pyfunction!(init_threads, m)?)?;
    m.add_function(wrap_pyfunction!(get_threads, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;

    // I/O functions
    m.add_function(wrap_pyfunction!(read_plink, m)?)?;
    m.add_function(wrap_pyfunction!(read_eigenstrat, m)?)?;
    m.add_function(wrap_pyfunction!(read_packedancestrymap, m)?)?;
    m.add_function(wrap_pyfunction!(read_vcf, m)?)?;
    m.add_function(wrap_pyfunction!(read_genotypes, m)?)?;
    m.add_function(wrap_pyfunction!(read_f2, m)?)?;
    m.add_function(wrap_pyfunction!(dataset_info, m)?)?;
    m.add_function(wrap_pyfunction!(list_populations, m)?)?;
    m.add_function(wrap_pyfunction!(suggest_populations, m)?)?;
    m.add_function(wrap_pyfunction!(convertf, m)?)?;
    m.add_function(wrap_pyfunction!(merge, m)?)?;
    m.add_function(wrap_pyfunction!(detect_format, m)?)?;

    // F-statistics functions
    m.add_function(wrap_pyfunction!(extract_f2, m)?)?;
    m.add_function(wrap_pyfunction!(f2_from_frequencies, m)?)?;
    geno_fstats::register(m)?;
    pca::register(m)?;
    kinship::register(m)?;
    sfs::register(m)?;
    admixture::register(m)?;
    qpadm_models::register(m)?;
    rcompat::register(m)?;
    rolloff::register(m)?;

    // qpAdm/qpWave functions
    m.add_function(wrap_pyfunction!(qpwave, m)?)?;
    m.add_function(wrap_pyfunction!(qpadm, m)?)?;
    m.add_function(wrap_pyfunction!(qpadm_multi, m)?)?;
    m.add_function(wrap_pyfunction!(qpf4ratio, m)?)?;
    m.add_function(wrap_pyfunction!(qpadm_rotate, m)?)?;
    m.add_function(wrap_pyfunction!(qpwave_pairs, m)?)?;

    // qpGraph functions
    m.add_function(wrap_pyfunction!(build_tree, m)?)?;
    m.add_function(wrap_pyfunction!(qpgraph, m)?)?;
    m.add_function(wrap_pyfunction!(find_graphs, m)?)?;
    m.add_function(wrap_pyfunction!(all_graphs, m)?)?;
    m.add_function(wrap_pyfunction!(graph_key, m)?)?;

    // DATES functions
    m.add_function(wrap_pyfunction!(dates, m)?)?;

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
    m.add_class::<PyRotationModel>()?;
    m.add_class::<PyQpAdmRotateResult>()?;
    m.add_class::<PyQpWavePairResult>()?;
    m.add_class::<PyQpWavePairsResult>()?;
    m.add_class::<PyAdmixGraph>()?;
    m.add_class::<PyScoredGraph>()?;
    m.add_class::<PyLdBin>()?;
    m.add_class::<PyDecayFit>()?;
    m.add_class::<PyDatesResult>()?;

    // Add version as module attribute
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
