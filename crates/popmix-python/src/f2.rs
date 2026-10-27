//! f2 blocks for Python (`F2Blocks`): f-statistics, FST and the block
//! arrays, read from and written to an f2 directory.

use crate::genotypes::PyAlleleFrequencies;
use crate::{compute, require_populations};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use popmix_core::fstats::{p_from_z, z_score};
use popmix_core::{BlockDefinition, F2Blocks};
use popmix_io::{read_f2_blocks, write_f2_blocks};
use pyo3::prelude::*;
use pyo3::types::PyDict;

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
        let refs: Vec<&str> = populations
            .iter()
            .map(std::string::String::as_str)
            .collect();
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
            .map(std::string::ToString::to_string)
            .collect()
    }

    /// f2(A, B).
    ///
    /// Returns:
    ///     dict with keys est, se, z, p, n_snps (as a row of a table), or
    ///     None if the pair has no data
    fn f2<'py>(
        &self,
        py: Python<'py>,
        pop_a: &str,
        pop_b: &str,
    ) -> PyResult<Option<Bound<'py, PyDict>>> {
        require_populations(&self.inner, &[pop_a, pop_b])?;
        self.inner
            .f2(pop_a, pop_b)
            .map(|(est, se, n_snps)| stat_dict(py, est, se, ("n_snps", n_snps)))
            .transpose()
    }

    /// f3(A; B, C) (admixture test): negative when A is admixed between B
    /// and C.
    ///
    /// Returns:
    ///     dict with keys est, se, z, p, n_blocks, or None
    fn f3<'py>(
        &self,
        py: Python<'py>,
        pop_a: &str,
        pop_b: &str,
        pop_c: &str,
    ) -> PyResult<Option<Bound<'py, PyDict>>> {
        require_populations(&self.inner, &[pop_a, pop_b, pop_c])?;
        self.inner
            .f3(pop_a, pop_b, pop_c)
            .map(|r| stat_dict(py, r.estimate, r.std_err, ("n_blocks", r.n_blocks)))
            .transpose()
    }

    /// f4(A, B; C, D): tests for gene flow between the (A, B) and (C, D)
    /// branches.
    ///
    /// Returns:
    ///     dict with keys est, se, z, p, n_blocks, or None
    fn f4<'py>(
        &self,
        py: Python<'py>,
        pop_a: &str,
        pop_b: &str,
        pop_c: &str,
        pop_d: &str,
    ) -> PyResult<Option<Bound<'py, PyDict>>> {
        require_populations(&self.inner, &[pop_a, pop_b, pop_c, pop_d])?;
        self.inner
            .f4(pop_a, pop_b, pop_c, pop_d)
            .map(|r| stat_dict(py, r.estimate, r.std_err, ("n_blocks", r.n_blocks)))
            .transpose()
    }

    /// FST (Hudson's fixation index) between two populations, as
    /// admixtools2's ``fst``: 0 for identical populations, 1 for fixed
    /// differences.
    ///
    /// Returns:
    ///     dict with keys est, se, n_snps, or None
    fn fst<'py>(
        &self,
        py: Python<'py>,
        pop_a: &str,
        pop_b: &str,
    ) -> PyResult<Option<Bound<'py, PyDict>>> {
        require_populations(&self.inner, &[pop_a, pop_b])?;
        self.inner
            .fst(pop_a, pop_b)
            .map(|r| {
                let d = PyDict::new(py);
                d.set_item("est", r.estimate)?;
                d.set_item("se", r.std_err)?;
                d.set_item("n_snps", r.n_snps)?;
                Ok(d)
            })
            .transpose()
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
    ) -> Option<Bound<'py, PyArray1<f64>>> {
        self.inner
            .get_block_f2_values(pop_a, pop_b)
            .map(|values| values.into_pyarray(py))
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
#[pyo3(signature = (af, *, block_size=None))]
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

/// A statistic as a dict with the keys of its table row, in order: est, se,
/// z and p (z NaN without an SE), then the count it rests on.
fn stat_dict<'py>(
    py: Python<'py>,
    est: f64,
    se: f64,
    count: (&str, usize),
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let z = z_score(est, se);
    d.set_item("est", est)?;
    d.set_item("se", se)?;
    d.set_item("z", z)?;
    d.set_item("p", p_from_z(z))?;
    d.set_item(count.0, count.1)?;
    Ok(d)
}
