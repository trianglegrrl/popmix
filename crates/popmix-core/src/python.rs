//! Errors as Python exceptions, for the bindings (feature `python`).
//!
//! Each error becomes the class of `popmix.errors` for its kind, so that a
//! binding can return a core error with `?`.

use pyo3::prelude::*;
use pyo3::types::PyType;

use crate::error::{ErrorKind, PopmixError};

/// The `popmix.errors` class an error is raised as.
pub fn exception_class(e: &PopmixError) -> &'static str {
    match (e, e.kind()) {
        (PopmixError::PopulationNotFound(_), _) => "PopulationNotFoundError",
        (_, ErrorKind::NotFound) => "MissingFileError",
        (_, ErrorKind::Io) => "FileError",
        (_, ErrorKind::Input) => "InputError",
        (_, ErrorKind::Computation) => "ComputationError",
        (_, ErrorKind::Internal) => "InternalError",
    }
}

/// A Python exception of `popmix.errors.<class>` with `message`.
pub fn popmix_exception(class: &str, message: String) -> PyErr {
    Python::attach(|py| {
        let class = py
            .import("popmix.errors")
            .and_then(|m| m.getattr(class))
            .and_then(|c| c.cast_into::<PyType>().map_err(PyErr::from));
        match class {
            Ok(class) => PyErr::from_type(class, message),
            Err(e) => e,
        }
    })
}

impl From<PopmixError> for PyErr {
    fn from(e: PopmixError) -> Self {
        popmix_exception(exception_class(&e), e.to_string())
    }
}

/// A `popmix.errors.InputError`: a bad argument, file or name.
pub fn input_error(message: impl Into<String>) -> PyErr {
    popmix_exception("InputError", message.into())
}

/// A `popmix.errors.MissingFileError`.
pub fn missing_file(message: impl Into<String>) -> PyErr {
    popmix_exception("MissingFileError", message.into())
}

/// A `popmix.errors.InternalError`: a bug in popmix.
pub fn internal_error(message: impl Into<String>) -> PyErr {
    popmix_exception("InternalError", message.into())
}
