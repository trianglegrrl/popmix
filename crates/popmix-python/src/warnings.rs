//! Warnings from the Rust core as Python warnings (UserWarning).

use popmix_core::{Result, Warning};
use popmix_io::VcfReader;

use pyo3::prelude::*;

/// Raise each warning as a Python UserWarning (an error, if the warnings
/// filter makes it one).
pub(crate) fn raise(warnings: Vec<Warning>) -> PyResult<()> {
    if warnings.is_empty() {
        return Ok(());
    }
    Python::attach(|py| {
        for w in warnings {
            PyErr::warn(
                py,
                &py.get_type::<pyo3::exceptions::PyUserWarning>(),
                &std::ffi::CString::new(w.to_string())?,
                1,
            )?;
        }
        Ok(())
    })
}

/// Run `f` with a sink for its warnings, then raise them.
pub(crate) fn with_warnings<T>(
    f: impl FnOnce(&mut dyn FnMut(Warning)) -> Result<T>,
) -> PyResult<T> {
    let mut warnings = Vec::new();
    let out = f(&mut |w| warnings.push(w));
    raise(warnings)?;
    Ok(out?)
}

/// Open a VCF, warning about the records skipped because they are not
/// biallelic SNPs: popmix's formats hold two single-base alleles.
pub(crate) fn open_vcf(path: impl AsRef<std::path::Path>) -> PyResult<VcfReader> {
    with_warnings(|warn| popmix_io::dataset::open_vcf(path, warn))
}
