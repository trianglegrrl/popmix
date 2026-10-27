//! Bindings for f4-, D- and f3-statistics computed from genotype files.

use numpy::{IntoPyArray, PyArray2};
use popmix_core::python::{internal_error, popmix_exception};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use popmix_core::{
    F3GenoOptions, F4Mode, GenoPopulations, GenoStats, GenotypeMatrix, GenotypeRows, MatrixRows,
    Sample, Snp,
};
use popmix_io::{PackedRows, SnpMajorRows};

use crate::genotypes::read_genotypes_of;
use crate::ploidy::{mapped, ploidy_rule, PloidyRule};

/// The SNPs, samples and row source for statistics on `populations`.
pub(crate) enum Source {
    Packed(PackedRows),
    Matrix(GenotypeMatrix, Vec<usize>, Vec<usize>),
}

/// Open `prefix` for the samples of `populations`: the row source (one row
/// per kept SNP, the selected samples in .ind order), the kept SNPs, the
/// samples' populations and ploidies, and the selected samples themselves.
pub(crate) fn open_source(
    prefix: &str,
    populations: &[String],
    auto_only: bool,
) -> PyResult<(Source, Vec<Snp>, GenoPopulations, Vec<Sample>)> {
    open_source_with(prefix, populations, auto_only, PloidyRule::Genome)
}

/// As `open_source`, with pseudo-haploid samples found by `rule`.
pub(crate) fn open_source_with(
    prefix: &str,
    populations: &[String],
    auto_only: bool,
    rule: PloidyRule,
) -> PyResult<(Source, Vec<Snp>, GenoPopulations, Vec<Sample>)> {
    open_source_filtered_with(
        prefix,
        populations,
        &|s: &Snp| !auto_only || s.is_autosomal(),
        rule,
    )
}

/// As `open_source`, keeping the SNPs for which `keep_snp` is true.
pub(crate) fn open_source_filtered(
    prefix: &str,
    populations: &[String],
    keep_snp: &dyn Fn(&Snp) -> bool,
) -> PyResult<(Source, Vec<Snp>, GenoPopulations, Vec<Sample>)> {
    open_source_filtered_with(prefix, populations, keep_snp, PloidyRule::Genome)
}

fn open_source_filtered_with(
    prefix: &str,
    populations: &[String],
    keep_snp: &dyn Fn(&Snp) -> bool,
    rule: PloidyRule,
) -> PyResult<(Source, Vec<Snp>, GenoPopulations, Vec<Sample>)> {
    let pop_refs: Vec<&str> = populations
        .iter()
        .map(std::string::String::as_str)
        .collect();
    let index = |pop: &str| populations.iter().position(|p| p == pop);
    let (samples, snps, file, matrix) = if let Some(mut m) = mapped(prefix)? {
        rule.apply_mapped(&m.file, &mut m.samples, &pop_refs);
        (m.samples, m.snps, Some(m.file), None)
    } else {
        let mut geno = read_genotypes_of(prefix, populations)?;
        rule.apply(&mut geno);
        (
            geno.samples().to_vec(),
            geno.snps().to_vec(),
            None,
            Some(geno),
        )
    };
    popmix_core::genotypes::check_sorted(&snps)?;

    let mut selected = Vec::new();
    let mut selected_samples = Vec::new();
    let mut sample_pop = Vec::new();
    let mut ploidy = Vec::new();
    for (i, s) in samples.iter().enumerate() {
        if let Some(p) = index(s.population.as_str()) {
            selected.push(i);
            selected_samples.push(s.clone());
            sample_pop.push(p);
            ploidy.push(s.ploidy.effective_ploidy());
        }
    }
    let missing: Vec<&str> = populations
        .iter()
        .filter(|p| !sample_pop.iter().any(|&k| populations[k] == **p))
        .map(std::string::String::as_str)
        .collect();
    if !missing.is_empty() {
        return Err(popmix_exception(
            "PopulationNotFoundError",
            format!("populations not in {prefix}: {}", missing.join(", ")),
        ));
    }

    let keep: Vec<bool> = snps.iter().map(keep_snp).collect();
    let kept: Vec<Snp> = snps
        .into_iter()
        .zip(&keep)
        .filter(|(_, &k)| k)
        .map(|(s, _)| s)
        .collect();
    let source = match (file, matrix) {
        (Some(file), _) => Source::Packed(SnpMajorRows::from_file(file, keep, selected)?),
        (None, Some(m)) => {
            let rows = (0..keep.len()).filter(|&i| keep[i]).collect();
            Source::Matrix(m, rows, selected)
        }
        (None, None) => unreachable!("a dataset is mapped or read"),
    };
    let pops = GenoPopulations {
        names: populations.to_vec(),
        sample_pop,
        ploidy,
    };
    Ok((source, kept, pops, selected_samples))
}

pub(crate) fn rows_of(source: &mut Source) -> Box<dyn GenotypeRows + '_> {
    match source {
        Source::Packed(rows) => Box::new(RowsRef(rows)),
        Source::Matrix(m, snps, samples) => {
            Box::new(MatrixRows::new(m, snps.clone(), samples.clone()))
        }
    }
}

struct RowsRef<'a>(&'a mut PackedRows);

impl GenotypeRows for RowsRef<'_> {
    fn next_row(&mut self, out: &mut [u8]) -> popmix_core::Result<()> {
        self.0.next_row(out)
    }

    fn next_rows(&mut self, out: &mut [u8], width: usize) -> popmix_core::Result<()> {
        self.0.next_rows(out, width)
    }
}

fn run<F>(source: Source, f: F) -> PyResult<GenoStats>
where
    F: FnOnce(&mut dyn GenotypeRows) -> popmix_core::Result<GenoStats>,
{
    match source {
        Source::Packed(mut rows) => f(&mut rows),
        Source::Matrix(m, snps, samples) => f(&mut MatrixRows::new(&m, snps, samples)),
    }
    .map_err(PyErr::from)
}

fn to_dict(py: Python<'_>, out: GenoStats) -> PyResult<Py<PyAny>> {
    let d = PyDict::new(py);
    let n_blocks = out.block_lengths.len();
    let col = |f: &dyn Fn(&popmix_core::GenoStat) -> f64| -> Vec<f64> {
        out.stats.iter().map(f).collect()
    };
    d.set_item("est", col(&|s| s.est).into_pyarray(py))?;
    d.set_item("se", col(&|s| s.se).into_pyarray(py))?;
    d.set_item("z", col(&|s| s.z).into_pyarray(py))?;
    d.set_item("p", col(&|s| s.p).into_pyarray(py))?;
    d.set_item(
        "n",
        out.stats
            .iter()
            .map(|s| s.n)
            .collect::<Vec<u64>>()
            .into_pyarray(py),
    )?;
    let blocks =
        |f: &dyn Fn(&popmix_core::GenoStat) -> Vec<f64>| -> PyResult<Bound<'_, PyArray2<f64>>> {
            let rows: Vec<Vec<f64>> = out.stats.iter().map(f).collect();
            let flat: Vec<f64> = rows.into_iter().flatten().collect();
            let arr = ndarray::Array2::from_shape_vec((out.stats.len(), n_blocks), flat)
                .map_err(|e| internal_error(e.to_string()))?;
            Ok(arr.into_pyarray(py))
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
#[pyo3(signature = (prefix, populations, quads, f4mode=true, block_size=0.05, auto_only=true, pseudohaploid="detect"))]
#[allow(clippy::too_many_arguments)]
pub fn f4_from_geno(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    quads: Vec<[usize; 4]>,
    f4mode: bool,
    block_size: f64,
    auto_only: bool,
    pseudohaploid: &str,
) -> PyResult<Py<PyAny>> {
    let rule = ploidy_rule(pseudohaploid)?;
    let (source, snps, pops, _) = open_source_with(prefix, &populations, auto_only, rule)?;
    warn_without_map(&snps, block_size)?;
    let mode = if f4mode { F4Mode::F4 } else { F4Mode::D };
    let out = crate::compute(py, || {
        run(source, |rows| {
            popmix_core::f4_from_genotypes(rows, &snps, &pops, &quads, mode, block_size)
        })
    })?;
    to_dict(py, out)
}

/// Warn when `snps` have no genetic map, so blocks fall back to 2 Mb (G5).
fn warn_without_map(snps: &[popmix_core::Snp], block_size: f64) -> PyResult<()> {
    let mut evidence = popmix_core::CoordinateEvidence::default();
    evidence.observe(snps);
    let blocks = popmix_core::BlockDefinition::new(block_size);
    crate::warnings::emit(|w| popmix_core::warnings::check_genetic_map(&evidence, &blocks, w))
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
#[pyo3(signature = (prefix, populations, triples, outgroupmode=false, apply_corr=true, block_size=0.05, auto_only=true, pseudohaploid="detect"))]
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
    pseudohaploid: &str,
) -> PyResult<Py<PyAny>> {
    let rule = ploidy_rule(pseudohaploid)?;
    let (source, snps, pops, _) = open_source_with(prefix, &populations, auto_only, rule)?;
    warn_without_map(&snps, block_size)?;
    let options = {
        let mut o = F3GenoOptions::default();
        o.apply_corr = apply_corr;
        o.outgroup_mode = outgroupmode;
        o
    };
    let out = crate::compute(py, || {
        run(source, |rows| {
            popmix_core::f3_from_genotypes(rows, &snps, &pops, &triples, options, block_size)
        })
    })?;
    to_dict(py, out)
}

/// Ploidy of each sample of ``populations``: 1 (pseudo-haploid) or 2, by
/// the ``pseudohaploid`` rule ("detect": no heterozygous genotype at any
/// SNP; "admixtools2": none among the first 1000 SNPs; "diploid").
#[pyfunction]
#[pyo3(signature = (prefix, populations, pseudohaploid="detect"))]
pub fn sample_ploidy(
    prefix: &str,
    populations: Vec<String>,
    pseudohaploid: &str,
) -> PyResult<Vec<(String, u8)>> {
    let rule = ploidy_rule(pseudohaploid)?;
    let pop_refs: Vec<&str> = populations
        .iter()
        .map(std::string::String::as_str)
        .collect();
    let samples = if let Some(mut m) = mapped(prefix)? {
        rule.apply_mapped(&m.file, &mut m.samples, &pop_refs);
        m.samples
    } else {
        let mut geno = read_genotypes_of(prefix, &populations)?;
        rule.apply(&mut geno);
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
