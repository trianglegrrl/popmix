//! Bindings for f4-, D- and f3-statistics computed from genotype files.

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use popmix_core::{
    F3GenoOptions, F4Mode, GenoPopulations, GenoStats, GenotypeMatrix, GenotypeRows, MatrixRows,
    Snp,
};
use popmix_io::{detect_ploidy_packedancestrymap, is_packed_format, PackedRows};

use crate::{read_genotypes_internal, to_py_err};

/// The SNPs, samples and row source for statistics on `populations`.
enum Source {
    Packed(PackedRows),
    Matrix(GenotypeMatrix, Vec<usize>, Vec<usize>),
}

fn open_source(
    prefix: &str,
    populations: &[String],
    auto_only: bool,
) -> PyResult<(Source, Vec<Snp>, GenoPopulations)> {
    let index = |pop: &str| populations.iter().position(|p| p == pop);
    let geno_path = format!("{prefix}.geno");
    let packed = std::path::Path::new(&geno_path).exists()
        && is_packed_format(&geno_path).map_err(to_py_err)?;
    let (samples, snps, matrix) = if packed {
        let mut samples =
            popmix_io::eigenstrat::read_ind_file(&format!("{prefix}.ind")).map_err(to_py_err)?;
        let snps =
            popmix_io::eigenstrat::read_snp_file(&format!("{prefix}.snp")).map_err(to_py_err)?;
        detect_ploidy_packedancestrymap(&geno_path, &mut samples, 1000, snps.len())
            .map_err(to_py_err)?;
        (samples, snps, None)
    } else {
        let mut geno = read_genotypes_internal(prefix)?;
        geno.detect_ploidy(1000);
        (geno.samples().to_vec(), geno.snps().to_vec(), Some(geno))
    };

    let mut selected = Vec::new();
    let mut sample_pop = Vec::new();
    let mut ploidy = Vec::new();
    for (i, s) in samples.iter().enumerate() {
        if let Some(p) = index(s.population.as_str()) {
            selected.push(i);
            sample_pop.push(p);
            ploidy.push(s.ploidy.effective_ploidy());
        }
    }
    let missing: Vec<&str> = populations
        .iter()
        .filter(|p| !sample_pop.iter().any(|&k| populations[k] == **p))
        .map(|p| p.as_str())
        .collect();
    if !missing.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "populations not in {prefix}: {}",
            missing.join(", ")
        )));
    }

    let keep: Vec<bool> = snps
        .iter()
        .map(|s| !auto_only || s.is_autosomal())
        .collect();
    let kept: Vec<Snp> = snps
        .iter()
        .zip(&keep)
        .filter(|(_, &k)| k)
        .map(|(s, _)| s.clone())
        .collect();
    let source = match matrix {
        None => Source::Packed(
            PackedRows::open(&geno_path, samples.len(), keep, selected).map_err(to_py_err)?,
        ),
        Some(m) => {
            let rows = (0..keep.len()).filter(|&i| keep[i]).collect();
            Source::Matrix(m, rows, selected)
        }
    };
    let pops = GenoPopulations {
        names: populations.to_vec(),
        sample_pop,
        ploidy,
    };
    Ok((source, kept, pops))
}

fn run<F>(source: Source, f: F) -> PyResult<GenoStats>
where
    F: FnOnce(&mut dyn GenotypeRows) -> popmix_core::Result<GenoStats>,
{
    match source {
        Source::Packed(mut rows) => f(&mut rows),
        Source::Matrix(m, snps, samples) => f(&mut MatrixRows::new(&m, snps, samples)),
    }
    .map_err(to_py_err)
}

fn to_dict(py: Python<'_>, out: GenoStats) -> PyResult<PyObject> {
    let d = PyDict::new_bound(py);
    let n_blocks = out.block_lengths.len();
    let col = |f: &dyn Fn(&popmix_core::GenoStat) -> f64| -> Vec<f64> {
        out.stats.iter().map(f).collect()
    };
    d.set_item("est", col(&|s| s.est).into_pyarray_bound(py))?;
    d.set_item("se", col(&|s| s.se).into_pyarray_bound(py))?;
    d.set_item("z", col(&|s| s.z).into_pyarray_bound(py))?;
    d.set_item("p", col(&|s| s.p).into_pyarray_bound(py))?;
    d.set_item(
        "n",
        out.stats
            .iter()
            .map(|s| s.n)
            .collect::<Vec<u64>>()
            .into_pyarray_bound(py),
    )?;
    let blocks =
        |f: &dyn Fn(&popmix_core::GenoStat) -> Vec<f64>| -> PyResult<Bound<'_, PyArray2<f64>>> {
            let rows: Vec<Vec<f64>> = out.stats.iter().map(f).collect();
            let flat: Vec<f64> = rows.into_iter().flatten().collect();
            let arr = ndarray::Array2::from_shape_vec((out.stats.len(), n_blocks), flat)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
            Ok(arr.into_pyarray_bound(py))
        };
    d.set_item("block_num", blocks(&|s| s.blocks.num.clone())?)?;
    d.set_item("block_den", blocks(&|s| s.blocks.den.clone())?)?;
    d.set_item(
        "block_n",
        blocks(&|s| s.blocks.n.iter().map(|&k| k as f64).collect())?,
    )?;
    d.set_item("block_lengths", out.block_lengths.clone())?;
    Ok(d.into())
}

/// f4- or D-statistics from a genotype prefix (admixtools2 ``qpdstat`` on
/// genotypes, ``allsnps = TRUE``).
///
/// Args:
///     prefix: Genotype prefix (PACKEDANCESTRYMAP is streamed; other
///         formats are read into memory)
///     populations: Population names that ``quads`` index into
///     quads: (pop1, pop2, pop3, pop4) index tuples
///     f4mode: True for f4, False for D-statistics
///     block_size: Jackknife block size in Morgans
///     auto_only: Use only autosomal SNPs
///
/// Returns:
///     dict of arrays: est, se, z, p, n (one per statistic), block_num,
///     block_den, block_n (statistics x blocks) and block_lengths
#[pyfunction]
#[pyo3(signature = (prefix, populations, quads, f4mode=true, block_size=0.05, auto_only=true))]
pub fn f4_from_geno(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    quads: Vec<[usize; 4]>,
    f4mode: bool,
    block_size: f64,
    auto_only: bool,
) -> PyResult<PyObject> {
    let (source, snps, pops) = open_source(prefix, &populations, auto_only)?;
    let mode = if f4mode { F4Mode::F4 } else { F4Mode::D };
    let out = py.allow_threads(|| {
        run(source, |rows| {
            popmix_core::f4_from_genotypes(rows, &snps, &pops, &quads, mode, block_size)
        })
    })?;
    to_dict(py, out)
}

/// f3-statistics from a genotype prefix (admixtools2 ``qp3pop`` on
/// genotypes, ``allsnps = TRUE``), pop1 the target.
///
/// Args:
///     prefix: Genotype prefix
///     populations: Population names that ``triples`` index into
///     triples: (pop1, pop2, pop3) index tuples
///     outgroupmode: Leave f3 unnormalized (otherwise it is divided by the
///         target's heterozygosity)
///     apply_corr: Subtract the target's sampling bias
///     block_size: Jackknife block size in Morgans
///     auto_only: Use only autosomal SNPs
///
/// Returns:
///     dict of arrays, as ``f4_from_geno``
#[pyfunction]
#[pyo3(signature = (prefix, populations, triples, outgroupmode=false, apply_corr=true, block_size=0.05, auto_only=true))]
#[allow(clippy::too_many_arguments)]
pub fn f3_from_geno(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    triples: Vec<[usize; 3]>,
    outgroupmode: bool,
    apply_corr: bool,
    block_size: f64,
    auto_only: bool,
) -> PyResult<PyObject> {
    let (source, snps, pops) = open_source(prefix, &populations, auto_only)?;
    let options = F3GenoOptions {
        apply_corr,
        outgroup_mode: outgroupmode,
    };
    let out = py.allow_threads(|| {
        run(source, |rows| {
            popmix_core::f3_from_genotypes(rows, &snps, &pops, &triples, options, block_size)
        })
    })?;
    to_dict(py, out)
}

/// Ploidy of each sample of ``populations``, as admixtools2 assigns it: 1
/// (pseudo-haploid) for a sample with no heterozygous genotype among the
/// first 1000 SNPs, otherwise 2.
#[pyfunction]
pub fn sample_ploidy(prefix: &str, populations: Vec<String>) -> PyResult<Vec<(String, u8)>> {
    let geno_path = format!("{prefix}.geno");
    let samples = if std::path::Path::new(&geno_path).exists()
        && is_packed_format(&geno_path).map_err(to_py_err)?
    {
        let mut samples =
            popmix_io::eigenstrat::read_ind_file(&format!("{prefix}.ind")).map_err(to_py_err)?;
        let n_snps = popmix_io::eigenstrat::read_snp_file(&format!("{prefix}.snp"))
            .map_err(to_py_err)?
            .len();
        detect_ploidy_packedancestrymap(&geno_path, &mut samples, 1000, n_snps)
            .map_err(to_py_err)?;
        samples
    } else {
        let mut geno = read_genotypes_internal(prefix)?;
        geno.detect_ploidy(1000);
        geno.samples().to_vec()
    };
    Ok(samples
        .into_iter()
        .filter(|s| populations.iter().any(|p| p == s.population.as_str()))
        .map(|s| (s.id.clone(), s.ploidy.effective_ploidy()))
        .collect())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(f4_from_geno, m)?)?;
    m.add_function(wrap_pyfunction!(f3_from_geno, m)?)?;
    m.add_function(wrap_pyfunction!(sample_ploidy, m)?)?;
    Ok(())
}
