//! Genotype data for Python: samples, SNPs, genotype matrices and allele
//! frequencies, and the readers of every genotype format.

use crate::warnings::{open_vcf, with_warnings};
use crate::{compute, require_sample_populations};
use crate::{geno_fstats, ploidy};
use numpy::{IntoPyArray, PyArray2};
use popmix_core::python::{input_error, missing_file};
use popmix_core::{AlleleFrequencies, GenotypeMatrix, Ploidy};
use popmix_io::tgeno::{is_tgeno_format, TgenoReader};
use popmix_io::{is_packed_format, EigenstratReader, PackedAncestryMapReader, PlinkReader};
use pyo3::prelude::*;

/// A sample (individual) with metadata.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PySample {
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
pub(crate) struct PySnp {
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
pub(crate) struct PyGenotypeMatrix {
    pub(crate) inner: GenotypeMatrix,
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
pub(crate) struct PyAlleleFrequencies {
    pub(crate) inner: AlleleFrequencies,
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
            inner: self.inner.filter_populations(&pop_refs)?,
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
pub(crate) fn read_plink(
    py: Python<'_>,
    prefix: &str,
    snp_range: Option<(usize, usize)>,
) -> PyResult<PyGenotypeMatrix> {
    compute(py, move || {
        let reader = PlinkReader::new(prefix)?;

        let geno = match snp_range {
            Some((start, end)) => reader.read_snp_range(start..end)?,
            None => reader.read_all()?,
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
pub(crate) fn read_eigenstrat(
    py: Python<'_>,
    prefix: &str,
    snp_range: Option<(usize, usize)>,
) -> PyResult<PyGenotypeMatrix> {
    compute(py, move || {
        let reader = EigenstratReader::new(prefix)?;

        let geno = match snp_range {
            Some((start, end)) => reader.read_snp_range(start..end)?,
            None => reader.read_all()?,
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
pub(crate) fn read_packedancestrymap(
    py: Python<'_>,
    prefix: &str,
    snp_range: Option<(usize, usize)>,
) -> PyResult<PyGenotypeMatrix> {
    compute(py, move || {
        let reader = PackedAncestryMapReader::new(prefix)?;

        let geno = match snp_range {
            Some((start, end)) => reader.read_snp_range(start..end)?,
            None => reader.read_all()?,
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
pub(crate) fn read_vcf(py: Python<'_>, path: &str) -> PyResult<PyGenotypeMatrix> {
    compute(py, move || {
        let reader = open_vcf(path)?;
        let geno = reader.read_all()?;
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
pub(crate) fn read_genotypes(
    py: Python<'_>,
    prefix: &str,
    snp_range: Option<(usize, usize)>,
    populations: Option<Vec<String>>,
) -> PyResult<PyGenotypeMatrix> {
    if let Some(pops) = populations {
        if snp_range.is_some() {
            return Err(input_error("give snp_range or populations, not both"));
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
pub(crate) fn read_genotypes_range(
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
        if is_tgeno_format(&geno_path)? {
            let reader = TgenoReader::new(prefix)?;
            let (start, end) = range;
            return Ok(PyGenotypeMatrix {
                inner: reader.read_snp_range(start..end)?,
            });
        }
        if is_packed_format(&geno_path)? {
            return read_packedancestrymap(py, prefix, snp_range);
        } else {
            return read_eigenstrat(py, prefix, snp_range);
        }
    }

    Err(missing_file(format!(
        "No genotype files found for prefix: {}",
        prefix
    )))
}

/// Every genotype of a dataset in any format, each sample pseudohaploid
/// when it has no heterozygous call (the default `pseudohaploid` rule).
/// Every SNP of the samples of `populations`, in file order: packed .geno
/// and PLINK .bed are read in place, other formats whole and then subset.
pub(crate) fn read_populations(prefix: &str, populations: &[String]) -> PyResult<GenotypeMatrix> {
    let (mut source, snps, _, samples) = geno_fstats::open_source(prefix, populations, false)?;
    let mut data = ndarray::Array2::<u8>::zeros((snps.len(), samples.len()));
    let mut rows = geno_fstats::rows_of(&mut source);
    for mut row in data.rows_mut() {
        rows.next_row(row.as_slice_mut().expect("row-major"))?;
    }
    drop(rows);
    GenotypeMatrix::new(data, samples, snps).map_err(PyErr::from)
}

pub(crate) fn read_genotypes_internal(prefix: &str) -> PyResult<GenotypeMatrix> {
    let mut geno = read_genotypes_any_format(prefix)?;
    ploidy::PloidyRule::Genome.apply(&mut geno);
    Ok(geno)
}

/// The genotypes of the samples of `populations` (or more): for TGENO only
/// their records are read, for other formats the whole file.
pub(crate) fn read_genotypes_of(prefix: &str, populations: &[String]) -> PyResult<GenotypeMatrix> {
    with_warnings(|warn| popmix_io::dataset::read_genotypes_of(prefix, populations, warn))
}

/// Every genotype of `prefix`, in any format.
pub(crate) fn read_genotypes_any_format(prefix: &str) -> PyResult<GenotypeMatrix> {
    with_warnings(|warn| popmix_io::dataset::read_genotypes(prefix, warn))
}

// ============================================================================
// DATES wrappers
// ============================================================================
