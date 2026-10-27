//! Bindings for admixtools2's qpAdm model-set functions.

use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::f2::PyF2Blocks;

/// admixtools2 ``qpadm_p``: one qpAdm model's fit at rank ``rnk``.
///
/// Returns a dict: f4rank, dof, chisq, p (None when dof is 0), feasible, and
/// weights (one per left population after the first, or None when the fit
/// does not determine them).
#[pyfunction(name = "qpadm_p")]
#[pyo3(signature = (f2_blocks, left, right, target=None, rnk=None, fudge=1e-4, constrained=false))]
#[allow(clippy::too_many_arguments)]
fn qpadm_p<'py>(
    py: Python<'py>,
    f2_blocks: &PyF2Blocks,
    left: Vec<String>,
    right: Vec<String>,
    target: Option<String>,
    rnk: Option<usize>,
    fudge: f64,
    constrained: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let left_refs: Vec<&str> = left.iter().map(|s| s.as_str()).collect();
    let right_refs: Vec<&str> = right.iter().map(|s| s.as_str()).collect();
    let r = popmix_core::qpadm_p(
        &f2_blocks.inner,
        &left_refs,
        &right_refs,
        target.as_deref(),
        rnk,
        fudge,
        constrained,
    )?;
    let d = PyDict::new(py);
    d.set_item("f4rank", r.f4rank)?;
    d.set_item("dof", r.dof)?;
    d.set_item("chisq", r.chi_sq)?;
    d.set_item("p", (!r.p_value.is_nan()).then_some(r.p_value))?;
    d.set_item("feasible", r.feasible)?;
    d.set_item("weights", r.weights)?;
    Ok(d)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(qpadm_p, m)?)?;
    Ok(())
}
