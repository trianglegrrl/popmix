//! Bindings for ALDER.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use popmix_core::alder::{AlderRun, CurveFit, Jackknifed, WeightedLdCurve};
use popmix_core::{AlderCompat, AlderOptions, PackedGenotypes};

use crate::geno_fstats::{open_source, rows_of};
use crate::to_py_err;

fn jk<'py>(py: Python<'py>, j: &Jackknifed) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("value", j.value)?;
    d.set_item("se", j.se())?;
    d.set_item("reps", j.reps.clone())?;
    Ok(d)
}

fn fit_dict<'py>(py: Python<'py>, f: &CurveFit) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("mindis", f.mindis)?;
    d.set_item("maxdis", f.maxdis)?;
    d.set_item("inter_chrom_affine", f.inter_chrom_affine)?;
    d.set_item("decay", jk(py, &f.decay)?)?;
    d.set_item("amp_tot", jk(py, &f.amp_tot)?)?;
    d.set_item("amp_exp", jk(py, &f.amp_exp)?)?;
    d.set_item("amp_aff", jk(py, &f.amp_aff)?)?;
    d.set_item("failures", f.failures)?;
    d.set_item("p", f.p_value(1.0))?;
    Ok(d)
}

fn curve_dict<'py>(py: Python<'py>, c: &WeightedLdCurve) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("d", c.d.clone())?;
    d.set_item("wld", c.wld.clone())?;
    d.set_item("count", c.count.clone())?;
    Ok(d)
}

fn run_dict<'py>(py: Python<'py>, r: &AlderRun) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("refs", r.refs.clone())?;
    d.set_item("fit_start", r.fit_start)?;
    d.set_item("curve", curve_dict(py, &r.curve)?)?;
    let jacks = PyList::empty_bound(py);
    for c in &r.jackknife_curves {
        jacks.append(curve_dict(py, c)?)?;
    }
    d.set_item("jackknife_curves", jacks)?;
    let fits = PyList::empty_bound(py);
    for f in &r.fits {
        fits.append(fit_dict(py, f)?)?;
    }
    d.set_item("fits", fits)?;
    d.set_item("fit_index", r.fit_index)?;
    Ok(d)
}

/// ALDER weighted LD with one or two reference populations.
///
/// Returns a dict: chromosomes, ld_corr_stops (Morgans, None with mindis),
/// run (the curve with ``refs`` weights: fit_start, curve, jackknife_curves,
/// fits, fit_index), oneref_runs (with two references), test, and
/// mixture_fraction (with one reference).
#[pyfunction(name = "alder")]
#[pyo3(signature = (prefix, test, refs, mindis=None, binsize=0.0005, maxdis=0.5, mincount=4, compat=None))]
#[allow(clippy::too_many_arguments)]
pub fn alder_py(
    py: Python<'_>,
    prefix: &str,
    test: String,
    refs: Vec<String>,
    mindis: Option<f64>,
    binsize: f64,
    maxdis: f64,
    mincount: usize,
    compat: Option<&str>,
) -> PyResult<PyObject> {
    let compat = match compat {
        None => AlderCompat::CORRECTED,
        Some("alder-1.03") => AlderCompat::V1_03,
        Some(other) => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "unknown compat {other:?} for alder; expected None or \"alder-1.03\""
            )))
        }
    };
    if refs.contains(&test) {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "the test population cannot also be a reference",
        ));
    }
    let populations: Vec<String> = std::iter::once(test.clone())
        .chain(refs.iter().cloned())
        .collect();
    let (mut source, snps, pops, samples) = open_source(prefix, &populations, true)?;
    let n = samples.len();
    let options = AlderOptions {
        binsize,
        maxdis,
        mincount,
        mindis,
        compat,
    };
    let out = py
        .allow_threads(|| {
            let mut rows = rows_of(&mut source);
            let all = PackedGenotypes::read(rows.as_mut(), n, snps.len())?;
            let of =
                |p: usize| -> Vec<usize> { (0..n).filter(|&s| pops.sample_pop[s] == p).collect() };
            let test_g = all.select_samples(&of(0));
            let ref_g: Vec<PackedGenotypes> = (1..populations.len())
                .map(|p| all.select_samples(&of(p)))
                .collect();
            let ref_refs: Vec<&PackedGenotypes> = ref_g.iter().collect();
            popmix_core::alder(&test_g, &ref_refs, &snps, &options)
        })
        .map_err(to_py_err)?;
    let d = PyDict::new_bound(py);
    d.set_item("chromosomes", out.chromosomes.clone())?;
    d.set_item("ld_corr_stops", out.ld_corr_stops.clone())?;
    d.set_item("run", run_dict(py, &out.run)?)?;
    let one = PyList::empty_bound(py);
    for r in &out.oneref_runs {
        one.append(run_dict(py, r)?)?;
    }
    d.set_item("oneref_runs", one)?;
    match &out.test {
        Some(t) => {
            let td = PyDict::new_bound(py);
            td.set_item("pretest", t.pretest)?;
            td.set_item("tworef_significant", t.tworef_significant)?;
            td.set_item("consistent", t.consistent)?;
            td.set_item("succeeds", t.succeeds)?;
            td.set_item("z", t.z)?;
            td.set_item("p", t.p)?;
            td.set_item(
                "decay_diffs",
                t.decay_diffs
                    .iter()
                    .map(|x| (x.z, x.percent))
                    .collect::<Vec<_>>(),
            )?;
            d.set_item("test", td)?;
        }
        None => d.set_item("test", py.None())?,
    }
    d.set_item("mixture_fraction", out.mixture_fraction)?;
    Ok(d.into())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(alder_py, m)?)?;
    Ok(())
}
