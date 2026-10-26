//! Bindings for ROLLOFF.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use popmix_core::{PackedGenotypes, RolloffBins, RolloffOptions, Snp};

use crate::geno_fstats::{open_source_filtered, rows_of};
use crate::to_py_err;

fn bins_dict<'py>(py: Python<'py>, b: &RolloffBins) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("d", b.d.clone())?;
    d.set_item("zscore", b.zscore.clone())?;
    d.set_item("regression", b.regression.clone())?;
    d.set_item("corr", b.corr.clone())?;
    d.set_item("count", b.count.clone())?;
    Ok(d)
}

/// ROLLOFF (ADMIXTOOLS 7.0.2) for ``admixed`` with references ``refs``
/// (two), then the fit of C + A exp(-m d) to the correlations from
/// ``lval`` to ``hval`` cM and its weighted jackknife over chromosomes (each
/// weighted by its SNPs in the .snp file), as expfit.sh.
///
/// Returns a dict: chromosomes, n_snps_used, bins, jackknife_bins (per
/// chromosome left out), fit (affine, amplitude, decay), jackknife_fits,
/// jackknife (mean, se), snp_weights.
#[pyfunction(name = "rolloff")]
#[pyo3(signature = (prefix, admixed, refs, binsize=0.0005, maxdis=0.05, runmode=1, zdipcorr=true,
                    chithresh=6.0, mincount=5, minparentcount=10, chromosomes=None, lval=0.5, hval=100.0))]
#[allow(clippy::too_many_arguments)]
pub fn rolloff_py(
    py: Python<'_>,
    prefix: &str,
    admixed: String,
    refs: Vec<String>,
    binsize: f64,
    maxdis: f64,
    runmode: u8,
    zdipcorr: bool,
    chithresh: f64,
    mincount: usize,
    minparentcount: usize,
    chromosomes: Option<Vec<String>>,
    lval: f64,
    hval: f64,
) -> PyResult<PyObject> {
    if refs.len() != 2 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "rolloff takes two reference populations, got {}",
            refs.len()
        )));
    }
    if refs.contains(&admixed) {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "the admixed population cannot also be a reference",
        ));
    }
    let chroms: Option<Vec<popmix_core::Chromosome>> = chromosomes.as_ref().map(|c| {
        c.iter()
            .map(|x| popmix_core::Chromosome::parse(x))
            .collect()
    });
    let keep =
        |s: &Snp| s.is_autosomal() && chroms.as_ref().is_none_or(|c| c.contains(&s.chromosome));
    let populations = vec![admixed.clone(), refs[0].clone(), refs[1].clone()];
    let (mut source, snps, pops, samples) = open_source_filtered(prefix, &populations, &keep)?;
    // jackknife weights: the SNPs of each chromosome in the .snp file
    let all_snps =
        popmix_io::eigenstrat::read_snp_file(&format!("{prefix}.snp")).map_err(to_py_err)?;
    let options = RolloffOptions {
        binsize,
        maxdis,
        runmode,
        zdipcorr,
        chithresh,
        mincount,
        minparentcount,
    };
    let n = samples.len();
    let out = crate::compute(py, || {
        let mut rows = rows_of(&mut source);
        let all = PackedGenotypes::read(rows.as_mut(), n, snps.len())?;
        let of = |p: usize| -> Vec<usize> { (0..n).filter(|&s| pops.sample_pop[s] == p).collect() };
        popmix_core::rolloff(
            &all.select_samples(&of(0)),
            &all.select_samples(&of(1)),
            &all.select_samples(&of(2)),
            &snps,
            &options,
        )
    })
    .map_err(to_py_err)?;
    let d_cm = |b: &RolloffBins| b.d.iter().map(|x| 100.0 * x).collect::<Vec<f64>>();
    let fit = popmix_core::fit_rolloff(&d_cm(&out.bins), &out.bins.corr, lval, hval);
    let jfits: Vec<popmix_core::RolloffFit> = out
        .jackknife
        .iter()
        .map(|b| popmix_core::fit_rolloff(&d_cm(b), &b.corr, lval, hval))
        .collect();
    let weight = |label: &str| {
        all_snps
            .iter()
            .filter(|s| s.chromosome.to_string() == label)
            .count() as f64
    };
    let blocks: Vec<(f64, f64)> = out
        .chromosomes
        .iter()
        .zip(&jfits)
        .map(|(c, f)| (weight(c), f.decay))
        .collect();
    let jack = popmix_core::dates::weighted_jackknife(fit.decay, &blocks).map_err(to_py_err)?;

    let d = PyDict::new_bound(py);
    d.set_item("chromosomes", out.chromosomes.clone())?;
    d.set_item("n_snps_used", out.n_snps_used)?;
    d.set_item("bins", bins_dict(py, &out.bins)?)?;
    let jl = PyList::empty_bound(py);
    for b in &out.jackknife {
        jl.append(bins_dict(py, b)?)?;
    }
    d.set_item("jackknife_bins", jl)?;
    d.set_item("fit", (fit.affine, fit.amplitude, fit.decay))?;
    d.set_item(
        "jackknife_fits",
        jfits
            .iter()
            .map(|f| (f.affine, f.amplitude, f.decay))
            .collect::<Vec<_>>(),
    )?;
    d.set_item("jackknife", jack)?;
    d.set_item(
        "snp_weights",
        blocks.iter().map(|b| b.0).collect::<Vec<_>>(),
    )?;
    Ok(d.into())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rolloff_py, m)?)?;
    Ok(())
}
