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

/// What `__reduce__` returns: the callable that rebuilds an object, and its
/// arguments.
type Reduce<'py, F, A> = (Bound<'py, F>, A);

/// A SNP's stored fields: ID, chromosome, position (bp, Morgans), alleles.
type SnpFields = (String, String, u32, f64, char, char);

/// A genotype matrix's parts: genotypes (SNPs x samples), samples, SNPs.
type MatrixParts<'py> = (Bound<'py, PyArray2<u8>>, Vec<PySample>, Vec<PySnp>);

/// A sample (individual): its ID, population and ploidy ("diploid" or
/// "pseudohaploid").
#[pyclass(name = "Sample", module = "popmix._core", frozen, eq, from_py_object)]
#[derive(Clone, PartialEq)]
pub(crate) struct PySample {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    population: String,
    #[pyo3(get)]
    ploidy: String,
}

#[pymethods]
impl PySample {
    #[new]
    #[pyo3(signature = (id, population, ploidy = "diploid"))]
    fn new(id: String, population: String, ploidy: &str) -> PyResult<Self> {
        if !matches!(ploidy, "diploid" | "pseudohaploid") {
            return Err(input_error(format!(
                "ploidy must be \"diploid\" or \"pseudohaploid\", got {ploidy:?}"
            )));
        }
        Ok(PySample {
            id,
            population,
            ploidy: ploidy.to_string(),
        })
    }

    fn __reduce__<'py>(
        slf: &Bound<'py, Self>,
    ) -> (Bound<'py, pyo3::types::PyType>, (String, String, String)) {
        let s = slf.get();
        (
            slf.get_type(),
            (s.id.clone(), s.population.clone(), s.ploidy.clone()),
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "Sample(id={:?}, population={:?}, ploidy={:?})",
            self.id, self.population, self.ploidy
        )
    }
}

impl PySample {
    pub(crate) fn to_core(&self) -> popmix_core::Sample {
        let mut s = popmix_core::Sample::new(&self.id, self.population.as_str());
        if self.ploidy == "pseudohaploid" {
            s.ploidy = Ploidy::Pseudohaploid;
        }
        s
    }
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

/// A SNP: its ID, chromosome, position (base pairs and centimorgans) and
/// alleles.
#[pyclass(name = "Snp", module = "popmix._core", frozen, eq, from_py_object)]
#[derive(Clone, PartialEq)]
pub(crate) struct PySnp {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    chromosome: String,
    #[pyo3(get)]
    position_bp: u32,
    /// kept in Morgans, as stored, so that pickling is exact
    position_morgans: f64,
    #[pyo3(get)]
    allele_ref: char,
    #[pyo3(get)]
    allele_alt: char,
}

#[pymethods]
impl PySnp {
    #[new]
    #[pyo3(signature = (id, chromosome, position_bp, position_cm = 0.0, allele_ref = 'N', allele_alt = 'N'))]
    fn new(
        id: String,
        chromosome: String,
        position_bp: u32,
        position_cm: f64,
        allele_ref: char,
        allele_alt: char,
    ) -> Self {
        PySnp {
            id,
            chromosome,
            position_bp,
            position_morgans: position_cm / 100.0,
            allele_ref,
            allele_alt,
        }
    }

    /// Genetic position in centimorgans.
    #[getter]
    fn position_cm(&self) -> f64 {
        self.position_morgans * 100.0
    }

    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<Reduce<'py, PyAny, SnpFields>> {
        Ok((
            py.import("popmix._core")?.getattr("_snp")?,
            (
                self.id.clone(),
                self.chromosome.clone(),
                self.position_bp,
                self.position_morgans,
                self.allele_ref,
                self.allele_alt,
            ),
        ))
    }

    fn __repr__(&self) -> String {
        format!(
            "Snp(id={:?}, chromosome={:?}, position_bp={}, position_cm={})",
            self.id,
            self.chromosome,
            self.position_bp,
            self.position_cm()
        )
    }
}

/// A SNP from its stored fields (genetic position in Morgans), for pickling.
#[pyfunction]
pub(crate) fn _snp(
    id: String,
    chromosome: String,
    position_bp: u32,
    position_morgans: f64,
    allele_ref: char,
    allele_alt: char,
) -> PySnp {
    PySnp {
        id,
        chromosome,
        position_bp,
        position_morgans,
        allele_ref,
        allele_alt,
    }
}

impl PySnp {
    pub(crate) fn to_core(&self) -> popmix_core::Snp {
        popmix_core::Snp::new(
            &self.id,
            popmix_core::Chromosome::parse(&self.chromosome),
            self.position_bp,
            self.position_morgans,
            self.allele_ref,
            self.allele_alt,
        )
    }
}

impl From<&popmix_core::Snp> for PySnp {
    fn from(s: &popmix_core::Snp) -> Self {
        PySnp {
            id: s.id.clone(),
            chromosome: s.chromosome.to_string(),
            position_bp: s.position_bp,
            position_morgans: s.position_morgans,
            allele_ref: s.allele_ref,
            allele_alt: s.allele_alt,
        }
    }
}

// ============================================================================
// GenotypeMatrix wrapper
// ============================================================================

/// A matrix of genotypes, SNPs by samples, with the samples and SNPs.
#[pyclass(name = "GenotypeMatrix", module = "popmix._core", frozen)]
pub(crate) struct PyGenotypeMatrix {
    pub(crate) inner: GenotypeMatrix,
}

#[pymethods]
impl PyGenotypeMatrix {
    /// A genotype matrix from its genotypes, SNPs by samples (0, 1, 2 copies
    /// of the counted allele, 3 missing), its samples and its SNPs.
    #[new]
    fn new(
        data: numpy::PyReadonlyArray2<'_, u8>,
        samples: Vec<PySample>,
        snps: Vec<PySnp>,
    ) -> PyResult<Self> {
        Ok(PyGenotypeMatrix {
            inner: GenotypeMatrix::new(
                data.as_array().to_owned(),
                samples.iter().map(PySample::to_core).collect(),
                snps.iter().map(PySnp::to_core).collect(),
            )?,
        })
    }

    fn __reduce__<'py>(
        slf: &Bound<'py, Self>,
    ) -> Reduce<'py, pyo3::types::PyType, MatrixParts<'py>> {
        let g = slf.get();
        let py = slf.py();
        (slf.get_type(), (g.to_numpy(py), g.samples(), g.snps()))
    }

    fn __repr__(&self) -> String {
        format!(
            "GenotypeMatrix({} SNPs x {} samples)",
            self.inner.n_snps(),
            self.inner.n_samples()
        )
    }

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
            .map(std::string::ToString::to_string)
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
                let pop_refs: Vec<&str> = pops.iter().map(std::string::String::as_str).collect();
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

/// Allele frequencies per population, SNPs by populations.
#[pyclass(name = "AlleleFrequencies", module = "popmix._core", frozen)]
pub(crate) struct PyAlleleFrequencies {
    pub(crate) inner: AlleleFrequencies,
}

/// Allele frequencies from their arrays, for pickling.
#[pyfunction]
pub(crate) fn _allele_frequencies(
    frequencies: numpy::PyReadonlyArray2<'_, f64>,
    counts: numpy::PyReadonlyArray2<'_, u32>,
    populations: Vec<String>,
    snps: Vec<PySnp>,
) -> PyResult<PyAlleleFrequencies> {
    Ok(PyAlleleFrequencies {
        inner: AlleleFrequencies::new(
            frequencies.as_array().to_owned(),
            counts.as_array().to_owned(),
            populations
                .iter()
                .map(|p| popmix_core::PopulationId::from(p.as_str()))
                .collect(),
            snps.iter().map(PySnp::to_core).collect(),
        )?,
    })
}

#[pymethods]
impl PyAlleleFrequencies {
    fn __reduce__<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<(Bound<'py, PyAny>, Bound<'py, pyo3::types::PyTuple>)> {
        let snps: Vec<PySnp> = self.inner.snps().iter().map(PySnp::from).collect();
        let args = pyo3::types::PyTuple::new(
            py,
            [
                self.inner.frequencies().clone().into_pyarray(py).into_any(),
                self.inner.counts().clone().into_pyarray(py).into_any(),
                pyo3::IntoPyObjectExt::into_bound_py_any(self.populations(), py)?,
                pyo3::IntoPyObjectExt::into_bound_py_any(snps, py)?,
            ],
        )?;
        Ok((
            py.import("popmix._core")?.getattr("_allele_frequencies")?,
            args,
        ))
    }

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
            .populations()
            .iter()
            .map(std::string::ToString::to_string)
            .collect()
    }

    /// Get allele frequencies as a numpy array.
    ///
    /// Returns:
    ///     2D numpy array of shape (n_snps, n_populations) with frequency values.
    ///     NaN indicates missing data.
    fn to_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.inner.frequencies().clone().into_pyarray(py)
    }

    /// Get observation counts as a numpy array.
    ///
    /// Returns:
    ///     2D numpy array of shape (n_snps, n_populations) with observation counts.
    fn counts_to_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<u32>> {
        self.inner.counts().clone().into_pyarray(py)
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
        let pop_refs: Vec<&str> = populations
            .iter()
            .map(std::string::String::as_str)
            .collect();
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
#[pyo3(signature = (prefix, *, snp_range=None))]
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
#[pyo3(signature = (prefix, *, snp_range=None))]
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
#[pyo3(signature = (prefix, *, snp_range=None))]
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
#[pyo3(signature = (prefix, *, snp_range=None, populations=None))]
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

    if let Some(path) = popmix_io::dataset::vcf_path(prefix) {
        // read whole: a VCF has no index to read a range by
        let inner = compute(py, || -> PyResult<GenotypeMatrix> {
            Ok(open_vcf(&path)?.read_all()?.select_snps(range.0..range.1)?)
        })?;
        return Ok(PyGenotypeMatrix { inner });
    }

    // Check for PLINK format (.bed file)
    if Path::new(&format!("{prefix}.bed")).exists() {
        return read_plink(py, prefix, snp_range);
    }

    // Check for EIGENSTRAT/PACKEDANCESTRYMAP format (.geno file)
    let geno_path = format!("{prefix}.geno");
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
        }
        return read_eigenstrat(py, prefix, snp_range);
    }

    Err(missing_file(format!(
        "No genotype files found for prefix: {prefix}"
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
