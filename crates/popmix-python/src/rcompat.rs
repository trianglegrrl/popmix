//! Bindings for the R-compatible numerical routines.

use pyo3::prelude::*;
use pyo3::types::PyDict;

use popmix_core::rcompat::{optim_lbfgsb, LbfgsbOptions, RRng};
use popmix_core::PopmixError;

/// R's ``set.seed(seed); runif(n)``.
#[pyfunction(name = "r_runif")]
fn r_runif(seed: i32, n: usize) -> Vec<f64> {
    RRng::set_seed(seed).runif(n)
}

/// ``reps`` draws of R's ``sample.int(n, size, replace = TRUE)`` after
/// ``set.seed(seed)``, 0-based: admixtools2's block bootstrap
/// (``boo_list``) draws them so.
#[pyfunction(name = "r_sample_replace")]
fn r_sample_replace(seed: i32, n: usize, size: usize, reps: usize) -> Vec<Vec<usize>> {
    let mut rng = RRng::set_seed(seed);
    (0..reps).map(|_| rng.sample_replace(n, size)).collect()
}

/// R's ``optim(par, fn, method = "L-BFGS-B", lower, upper, control)``
/// with finite-difference gradients.
///
/// Returns a dict: par, value, evaluations, convergence, message.
#[pyfunction(name = "r_optim_lbfgsb")]
#[pyo3(signature = (fn_, par, lower, upper, fnscale=1.0, maxit=100, lmm=5, factr=1e7, pgtol=0.0, ndeps=1e-3))]
#[allow(clippy::too_many_arguments)]
fn r_optim_lbfgsb<'py>(
    py: Python<'py>,
    fn_: Py<PyAny>,
    par: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    fnscale: f64,
    maxit: usize,
    lmm: usize,
    factr: f64,
    pgtol: f64,
    ndeps: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let mut f = |x: &[f64]| -> popmix_core::Result<f64> {
        fn_.call1(py, (x.to_vec(),))
            .and_then(|v| v.extract::<f64>(py))
            .map_err(|e| PopmixError::numerical("optim", format!("the objective failed: {e}")))
    };
    let opts = {
        let mut o = LbfgsbOptions::default();
        o.lmm = lmm;
        o.factr = factr;
        o.pgtol = pgtol;
        o.maxit = maxit;
        o.fnscale = fnscale;
        o.ndeps = ndeps;
        o
    };
    let r = optim_lbfgsb(&mut f, None, &par, &lower, &upper, &opts)?;
    let d = PyDict::new(py);
    d.set_item("par", r.par)?;
    d.set_item("value", r.value)?;
    d.set_item("evaluations", r.evaluations)?;
    d.set_item("convergence", r.convergence)?;
    d.set_item("message", r.message)?;
    Ok(d)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(r_runif, m)?)?;
    m.add_function(wrap_pyfunction!(r_sample_replace, m)?)?;
    m.add_function(wrap_pyfunction!(r_optim_lbfgsb, m)?)?;
    Ok(())
}
