//! f2 blocks for Python (`F2Blocks`): f-statistics, FST and the block
//! arrays, read from and written to an f2 directory.

use crate::genotypes::PyAlleleFrequencies;
use crate::{compute, require_populations};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use popmix_core::fstats::{p_from_z, z_score};
use popmix_core::{BlockDefinition, F2Blocks};
use popmix_io::{read_f2_blocks, write_f2_blocks};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use pyo3::IntoPyObjectExt;

/// Load f2 blocks saved by ``F2Blocks.save``.
///
/// Directories of every storage version are read; those before version 4
/// hold no FST blocks, so ``fst`` is approximated from the f2 blocks.
///
/// Parameters
/// ----------
/// path
///     The directory, holding ``metadata.json`` and the ``block_*.npy``
///     arrays.
///
/// Returns
/// -------
/// F2Blocks
///     The f2 blocks, ready for f-statistics.
#[pyfunction]
pub(crate) fn read_f2(py: Python<'_>, path: &str) -> PyResult<PyF2Blocks> {
    let inner = compute(py, move || -> PyResult<F2Blocks> {
        Ok(read_f2_blocks(std::path::Path::new(path))?)
    })?;
    crate::warnings::checked_blocks(inner)
}

/// Per-block f2 statistics for every pair of populations.
///
/// The data structure f-statistics are computed from, with block jackknife
/// standard errors. Made by ``extract_f2``, ``f2_from_frequencies``,
/// ``read_f2`` or ``F2Blocks.from_block_means``.
#[pyclass(name = "F2Blocks", module = "popmix._core", frozen)]
pub(crate) struct PyF2Blocks {
    pub(crate) inner: F2Blocks,
}

impl PyF2Blocks {
    /// Warn for each of `pops` with at most one allele per SNP, whose
    /// sampling bias `statistic` cannot correct (TASKS.md G7).
    fn warn_one_allele(&self, pops: &[&str], statistic: &str) -> PyResult<()> {
        crate::warnings::emit(|w| {
            for pop in pops.iter().filter(|p| self.inner.has_one_allele(p)) {
                w(popmix_core::Warning::OneAllele {
                    population: (*pop).to_string(),
                    statistic: statistic.to_string(),
                });
            }
        })
    }
}

#[pymethods]
impl PyF2Blocks {
    fn __reduce__<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyTuple>)> {
        let f = &self.inner;
        let array = |a: &ndarray::Array2<f64>| a.clone().into_pyarray(py).into_any();
        let fst = match f.fst_blocks() {
            Some(fst) => (
                fst.block_lengths().to_vec(),
                array(fst.num()),
                array(fst.den()),
            )
                .into_bound_py_any(py)?,
            None => py.None().into_bound(py),
        };
        let args = PyTuple::new(
            py,
            [
                array(f.block_sums()),
                array(f.block_raw_sums()),
                array(f.block_weights()),
                array(f.block_het_sums()),
                array(f.block_het_weights()),
                self.populations().into_bound_py_any(py)?,
                f.block_lengths().to_vec().into_bound_py_any(py)?,
                fst,
            ],
        )?;
        Ok((py.import("popmix._core")?.getattr("_f2_blocks")?, args))
    }

    fn __repr__(&self) -> String {
        format!(
            "F2Blocks({} populations, {} blocks)",
            self.inner.populations().len(),
            self.inner.n_blocks()
        )
    }

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

    /// Keep only some populations (admixtools2 ``extract_f2_subset``).
    ///
    /// Parameters
    /// ----------
    /// populations
    ///     The populations to keep, in the order wanted. An unknown one
    ///     raises ``PopulationNotFoundError``.
    ///
    /// Returns
    /// -------
    /// F2Blocks
    ///     New f2 blocks of those populations, FST blocks kept.
    fn subset(&self, populations: Vec<String>) -> PyResult<PyF2Blocks> {
        let refs: Vec<&str> = populations
            .iter()
            .map(std::string::String::as_str)
            .collect();
        Ok(PyF2Blocks {
            inner: self.inner.subset_populations(&refs)?,
        })
    }

    /// SNPs with data for both populations in each block.
    ///
    /// These are the counts admixtools2 stores with its f2 blocks.
    ///
    /// Parameters
    /// ----------
    /// pop_a, pop_b
    ///     The two populations.
    ///
    /// Returns
    /// -------
    /// list of float or None
    ///     One count per block, or None if a population is unknown.
    fn counts_per_block(&self, pop_a: &str, pop_b: &str) -> Option<Vec<f64>> {
        self.inner.block_counts(pop_a, pop_b)
    }

    /// F2 blocks from each pair's per-block f2 mean and SNP count.
    ///
    /// This is how admixtools2 stores them. FST is not available from
    /// these (its sums are NaN).
    ///
    /// Parameters
    /// ----------
    /// populations
    ///     The population names, in index order.
    /// block_lengths
    ///     Number of SNPs in each block (the jackknife weights).
    /// means
    ///     ``means[i][j]``: the per-block f2 means of populations i and j
    ///     (only i <= j is read).
    /// counts
    ///     ``counts[i][j]``: the per-block SNP counts of populations i and j,
    ///     likewise. Shapes that do not match raise ``InputError``.
    ///
    /// Returns
    /// -------
    /// F2Blocks
    ///     The f2 blocks.
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
        crate::warnings::checked_blocks(F2Blocks::from_block_means(
            pops,
            block_lengths,
            &means,
            &counts,
        )?)
    }

    /// Select blocks, as admixtools2 ``f2_blocks[, , blocks]``.
    ///
    /// For resampling such as the bootstrap. FST blocks are not carried
    /// over.
    ///
    /// Parameters
    /// ----------
    /// blocks
    ///     Block indices from 0, in the order wanted; repeats allowed. An
    ///     empty list or an index out of range raises ``InputError``.
    ///
    /// Returns
    /// -------
    /// F2Blocks
    ///     New f2 blocks with those blocks.
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

    /// Per-block FST of a pair of populations, on the FST blocks.
    ///
    /// Parameters
    /// ----------
    /// pop_a, pop_b
    ///     The two populations.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray or None
    ///     float64 array of shape (len(fst_block_lengths),), NaN for a block
    ///     without data; None if a population is unknown or the blocks carry
    ///     no FST data.
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

    /// f2(A, B), with a weighted block jackknife standard error.
    ///
    /// f2(A, A) is 0. An unknown population raises
    /// ``PopulationNotFoundError``.
    ///
    /// Parameters
    /// ----------
    /// pop_a, pop_b
    ///     The two populations.
    ///
    /// Returns
    /// -------
    /// dict or None
    ///     ``est``, ``se``, ``z``, ``p`` (z is NaN without an SE) and
    ///     ``n_snps`` (a row of a table), or None if no block has data.
    fn f2<'py>(
        &self,
        py: Python<'py>,
        pop_a: &str,
        pop_b: &str,
    ) -> PyResult<Option<Bound<'py, PyDict>>> {
        require_populations(&self.inner, &[pop_a, pop_b])?;
        if pop_a != pop_b {
            self.warn_one_allele(&[pop_a, pop_b], "f2 with it")?;
        }
        self.inner
            .f2(pop_a, pop_b)
            .map(|(est, se, n_snps)| stat_dict(py, est, se, ("n_snps", n_snps)))
            .transpose()
    }

    /// f3(A; B, C), the admixture test: negative when A is admixed
    /// between B and C.
    ///
    /// An unknown population raises ``PopulationNotFoundError``.
    ///
    /// Parameters
    /// ----------
    /// pop_a
    ///     The target A.
    /// pop_b, pop_c
    ///     The sources B and C.
    ///
    /// Returns
    /// -------
    /// dict or None
    ///     ``est``, ``se``, ``z``, ``p`` (z is NaN without an SE) and
    ///     ``n_blocks`` (blocks with data), or None if no block has data.
    fn f3<'py>(
        &self,
        py: Python<'py>,
        pop_a: &str,
        pop_b: &str,
        pop_c: &str,
    ) -> PyResult<Option<Bound<'py, PyDict>>> {
        require_populations(&self.inner, &[pop_a, pop_b, pop_c])?;
        self.warn_one_allele(&[pop_a], "f3 with it as the target")?;
        self.inner
            .f3(pop_a, pop_b, pop_c)
            .map(|r| stat_dict(py, r.estimate, r.std_err, ("n_blocks", r.n_blocks)))
            .transpose()
    }

    /// f4(A, B; C, D): tests for gene flow between the (A, B) and (C, D)
    /// branches.
    ///
    /// An unknown population raises ``PopulationNotFoundError``.
    ///
    /// Parameters
    /// ----------
    /// pop_a, pop_b, pop_c, pop_d
    ///     The populations A, B, C and D.
    ///
    /// Returns
    /// -------
    /// dict or None
    ///     ``est``, ``se``, ``z``, ``p`` (z is NaN without an SE) and
    ///     ``n_blocks`` (blocks with data), or None if no block has data.
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

    /// FST (Hudson's fixation index) between two populations.
    ///
    /// As admixtools2's ``fst``: 0 for identical populations, 1 for fixed
    /// differences, from the FST blocks; f2 blocks saved before storage
    /// version 4 have none, and FST is then approximated from the f2
    /// blocks. An unknown population raises ``PopulationNotFoundError``.
    ///
    /// Parameters
    /// ----------
    /// pop_a, pop_b
    ///     The two populations.
    ///
    /// Returns
    /// -------
    /// dict or None
    ///     ``est``, ``se`` and ``n_snps``, or None if no block has data.
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

    /// Get the f2 estimates of every pair as a matrix.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     Symmetric float64 array of shape (n_populations, n_populations),
    ///     in ``populations`` order: the ``f2`` estimates, 0 on the diagonal,
    ///     NaN where a pair has no SNPs.
    fn f2_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.inner.f2_matrix().into_pyarray(py)
    }

    /// Get the per-block f2 means of a pair of populations.
    ///
    /// Parameters
    /// ----------
    /// pop_a, pop_b
    ///     The two populations.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray or None
    ///     float64 array of shape (n_blocks,), NaN for a block without
    ///     SNPs and all 0 for a population with itself; None if a population
    ///     is unknown.
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

    /// Save the f2 blocks to a directory, for ``read_f2``.
    ///
    /// Writes ``metadata.json`` (populations, block lengths, storage
    /// version) and the ``block_*.npy`` arrays: f2 and raw sums, weights,
    /// heterozygosity sums and weights, and the FST blocks when present.
    ///
    /// Parameters
    /// ----------
    /// path
    ///     The directory, created if needed.
    /// source
    ///     Path of the genotype data, recorded in the metadata, or None (the
    ///     default).
    ///
    /// Returns
    /// -------
    /// None
    #[pyo3(signature = (path, source=None))]
    fn save(&self, path: &str, source: Option<&str>) -> PyResult<()> {
        write_f2_blocks(&self.inner, std::path::Path::new(path), source).map_err(PyErr::from)
    }
}

// ============================================================================
// F-statistics functions
// ============================================================================

/// Compute f2 blocks from allele frequencies.
///
/// SNPs are chosen as ``extract_f2`` does by default (no missing data,
/// polymorphic only).
///
/// Parameters
/// ----------
/// af
///     The allele frequencies (``AlleleFrequencies``).
/// block_size
///     Jackknife block size in Morgans (default 0.05, 5 cM).
///
/// Returns
/// -------
/// F2Blocks
///     The f2 blocks of the populations of ``af``. FST is approximated
///     from these, as for f2 blocks without FST blocks.
#[pyfunction]
#[pyo3(signature = (af, *, block_size=0.05))]
pub(crate) fn f2_from_frequencies(
    py: Python<'_>,
    af: &PyAlleleFrequencies,
    block_size: f64,
) -> PyResult<PyF2Blocks> {
    let mut evidence = popmix_core::CoordinateEvidence::default();
    evidence.observe(af.inner.snps());
    crate::warnings::emit(|w| {
        popmix_core::warnings::check_genetic_map(&evidence, &BlockDefinition::new(block_size), w);
    })?;
    let inner = compute(py, move || {
        let block_def = BlockDefinition::new(block_size);
        F2Blocks::from_allele_frequencies(&af.inner, &block_def)
    });
    crate::warnings::checked_blocks(inner)
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

type Arr<'a> = numpy::PyReadonlyArray2<'a, f64>;

/// F2 blocks from their arrays, for pickling (`F2Blocks.__reduce__`).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub(crate) fn _f2_blocks(
    block_sums: Arr<'_>,
    block_raw_sums: Arr<'_>,
    block_weights: Arr<'_>,
    block_het_sums: Arr<'_>,
    block_het_weights: Arr<'_>,
    populations: Vec<String>,
    block_lengths: Vec<usize>,
    fst: Option<(Vec<usize>, Arr<'_>, Arr<'_>)>,
) -> PyResult<PyF2Blocks> {
    let fst = match fst {
        Some((lengths, num, den)) => Some(popmix_core::fstats::FstBlocks::new(
            lengths,
            num.as_array().to_owned(),
            den.as_array().to_owned(),
        )?),
        None => None,
    };
    let inner = F2Blocks::from_arrays(popmix_core::F2Arrays {
        block_sums: block_sums.as_array().to_owned(),
        block_raw_sums: block_raw_sums.as_array().to_owned(),
        block_weights: block_weights.as_array().to_owned(),
        block_het_sums: block_het_sums.as_array().to_owned(),
        block_het_weights: block_het_weights.as_array().to_owned(),
        populations: populations
            .iter()
            .map(|p| popmix_core::PopulationId::from(p.as_str()))
            .collect(),
        block_lengths,
        fst,
    })?;
    Ok(PyF2Blocks { inner })
}
