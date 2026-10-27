//! The `pseudohaploid` option and the datasets read in place, for the
//! bindings: `popmix_io::dataset`, with its warnings raised in Python.

use pyo3::prelude::*;

pub(crate) use popmix_io::dataset::{Mapped, PloidyRule};

use crate::warnings::with_warnings;

/// `prefix` as a dataset read in place, or None if its format cannot be.
pub(crate) fn mapped(prefix: &str) -> PyResult<Option<Mapped>> {
    with_warnings(|warn| popmix_io::dataset::open_mapped(prefix, warn))
}

/// Parse the `pseudohaploid` option.
pub(crate) fn ploidy_rule(name: &str) -> PyResult<PloidyRule> {
    Ok(name.parse()?)
}
