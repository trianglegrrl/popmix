//! Bindings for ROLLOFF.

use popmix_core::python::input_error;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use popmix_core::{PackedGenotypes, RolloffBins, RolloffOptions, Snp};

use crate::geno_fstats::{open_source_filtered, rows_of};

fn bins_dict<'py>(py: Python<'py>, b: &RolloffBins) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
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
) -> PyResult<Py<PyAny>> {
    if refs.len() != 2 {
        return Err(input_error(format!(
            "rolloff takes two reference populations, got {}",
            refs.len()
        )));
    }
    if refs.contains(&admixed) {
        return Err(input_error(
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
    let all_snps = popmix_io::eigenstrat::read_snp_file(format!("{prefix}.snp"))?;
    let options = {
        let mut o = RolloffOptions::default();
        o.binsize = binsize;
        o.maxdis = maxdis;
        o.runmode = runmode;
        o.zdipcorr = zdipcorr;
        o.chithresh = chithresh;
        o.mincount = mincount;
        o.minparentcount = minparentcount;
        o
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
    })?;
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
    let jack = popmix_core::dates::weighted_jackknife(fit.decay, &blocks)?;
    // the two-pulse check (TASKS.md H8), on the bins of the fitting range
    let in_range = |b: &RolloffBins| -> Vec<f64> {
        b.d.iter()
            .zip(&b.corr)
            .filter(|(d, _)| (lval..=hval).contains(&(100.0 * **d)))
            .map(|(_, y)| *y)
            .collect()
    };
    let two = |b: &RolloffBins| popmix_core::dates::fit_two_pulse_values(&in_range(b), binsize);
    let pulse_check = match (
        two(&out.bins),
        out.jackknife.iter().map(two).collect::<Option<Vec<_>>>(),
    ) {
        (Some(full), Some(left_out)) if left_out.len() >= 2 => {
            let weighted: Vec<_> = blocks.iter().map(|b| b.0).zip(left_out).collect();
            Some(popmix_core::dates::check_pulses(&full, &weighted)?)
        }
        _ => None,
    };
    crate::warnings::emit(|w| {
        popmix_core::warnings::check_chromosomes("ROLLOFF", out.chromosomes.len(), w);
    })?;

    let d = PyDict::new(py);
    d.set_item("chromosomes", out.chromosomes.clone())?;
    d.set_item("n_snps_used", out.n_snps_used)?;
    d.set_item("bins", bins_dict(py, &out.bins)?)?;
    let jl = PyList::empty(py);
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
    match pulse_check {
        Some(c) => {
            let p = PyDict::new(py);
            p.set_item("generations", (c.fit.generations[0], c.fit.generations[1]))?;
            p.set_item("generations_se", (c.generations_se[0], c.generations_se[1]))?;
            p.set_item("amplitudes", (c.fit.amplitudes[0], c.fit.amplitudes[1]))?;
            p.set_item("amplitudes_se", (c.amplitudes_se[0], c.amplitudes_se[1]))?;
            p.set_item("difference_se", c.difference_se)?;
            p.set_item("rss", c.fit.rss)?;
            p.set_item("multiple_pulses", c.multiple_pulses)?;
            d.set_item("pulse_check", p)?;
        }
        None => d.set_item("pulse_check", py.None())?,
    }
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
