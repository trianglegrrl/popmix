//! f2 blocks for Python (`F2Blocks`): f-statistics, FST and the block
//! arrays, read from and written to an f2 directory.

use crate::genotypes::PyAlleleFrequencies;
use crate::{compute, require_populations};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use popmix_core::{BlockDefinition, F2Blocks};
use popmix_io::{read_f2_blocks, write_f2_blocks};
use pyo3::prelude::*;
use pyo3::IntoPyObjectExt;
use std::collections::{BTreeMap, HashMap};

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
pub(crate) fn read_f2(py: Python<'_>, path: &str) -> PyResult<PyF2Blocks> {
    compute(py, move || {
        let inner = read_f2_blocks(std::path::Path::new(path))?;
        Ok(PyF2Blocks { inner })
    })
}

/// Precomputed F2 statistics for all population pairs.
///
/// This is the main data structure for computing f-statistics.
/// Create using `extract_f2()` or from allele frequencies.
#[pyclass]
pub(crate) struct PyF2Blocks {
    pub(crate) inner: F2Blocks,
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
        self.inner.n_blocks()
    }

    /// These populations only, in this order (admixtools2
    /// ``extract_f2_subset``); FST blocks are kept.
    fn subset(&self, populations: Vec<String>) -> PyResult<PyF2Blocks> {
        let refs: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
        Ok(PyF2Blocks {
            inner: self.inner.subset_populations(&refs)?,
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
            inner: F2Blocks::from_block_means(pops, block_lengths, &means, &counts)?,
        })
    }

    /// These blocks, in this order and with repeats, as admixtools2
    /// ``f2_blocks[, , blocks]`` (0-based). FST blocks are not carried over.
    fn select_blocks(&self, blocks: Vec<usize>) -> PyResult<PyF2Blocks> {
        Ok(PyF2Blocks {
            inner: self.inner.select_blocks(&blocks)?,
        })
    }

    /// Number of SNPs in each block (the jackknife weights).
    #[getter]
    fn block_lengths(&self) -> Vec<usize> {
        self.inner.block_lengths().to_vec()
    }

    /// Number of SNPs in each FST block (FST uses its own SNP set and blocks,
    /// as admixtools2). Empty if the blocks carry no FST data (f2 directories
    /// saved before storage version 4).
    #[getter]
    fn fst_block_lengths(&self) -> Vec<usize> {
        self.inner
            .fst_blocks()
            .map(|f| f.block_lengths().to_vec())
            .unwrap_or_default()
    }

    /// Per-block FST for a pair of populations, or None if unavailable.
    fn fst_per_block<'py>(
        &self,
        py: Python<'py>,
        pop_a: &str,
        pop_b: &str,
    ) -> Option<Bound<'py, PyArray1<f64>>> {
        let fst = self.inner.fst_blocks()?;
        let pair = self
            .inner
            .pair_index(self.inner.pop_index(pop_a)?, self.inner.pop_index(pop_b)?);
        Some(ndarray::Array1::from(fst.block_values(pair)).into_pyarray(py))
    }

    /// Get population names.
    #[getter]
    fn populations(&self) -> Vec<String> {
        self.inner
            .populations()
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
        write_f2_blocks(&self.inner, std::path::Path::new(path), source).map_err(PyErr::from)
    }
}

// ============================================================================
// F-statistics functions
// ============================================================================

/// Compute F2 blocks from allele frequencies.
#[pyfunction]
#[pyo3(signature = (af, block_size=None))]
pub(crate) fn f2_from_frequencies(
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
