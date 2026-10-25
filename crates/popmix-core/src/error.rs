//! Error types for popmix.
//!
//! All errors are explicit and propagated - NO exception swallowing.

use thiserror::Error;

/// Result type alias for popmix operations.
pub type Result<T> = std::result::Result<T, PopmixError>;

/// Comprehensive error type for all popmix operations.
///
/// Errors are categorized by domain and carry sufficient context for debugging.
#[derive(Error, Debug)]
pub enum PopmixError {
    // === I/O Errors ===
    #[error("File not found: {path}")]
    FileNotFound { path: String },

    #[error("Invalid file format in '{file}': {details}")]
    InvalidFileFormat { file: String, details: String },

    #[error("Unexpected file size for '{file}': expected {expected} bytes, got {actual}")]
    UnexpectedFileSize {
        file: String,
        expected: u64,
        actual: u64,
    },

    #[error("Failed to parse '{file}' at line {line}: {details}")]
    ParseError {
        file: String,
        line: usize,
        details: String,
    },

    // === Data Validation Errors ===
    #[error("Population not found: '{0}'")]
    PopulationNotFound(String),

    #[error("Duplicate population: '{0}'")]
    DuplicatePopulation(String),

    #[error("Insufficient populations: need at least {required}, got {provided}")]
    InsufficientPopulations { required: usize, provided: usize },

    #[error("No SNPs remaining after filtering")]
    NoSnpsRemaining,

    #[error("Insufficient data for {context}: need at least {required}, got {actual}")]
    InsufficientData {
        context: String,
        required: usize,
        actual: usize,
    },

    #[error("Dimension mismatch in {context}: expected {expected}, got {actual}")]
    DimensionMismatch {
        context: String,
        expected: String,
        actual: String,
    },

    #[error("Invalid population subset: {details}")]
    InvalidSubset { details: String },

    // === Numerical Errors ===
    #[error("Matrix is singular or near-singular in {context}")]
    SingularMatrix { context: String },

    #[error("Optimization failed to converge after {iterations} iterations in {context}")]
    ConvergenceFailure { iterations: usize, context: String },

    #[error("No feasible solution in {context}")]
    Infeasible { context: String },

    #[error("Numerical overflow in {context}")]
    NumericalOverflow { context: String },

    #[error("Invalid value: {context} - {details}")]
    InvalidValue { context: String, details: String },

    // === Graph Errors ===
    #[error("Invalid graph structure: {details}")]
    InvalidGraph { details: String },

    #[error("Graph has no root node")]
    NoRootNode,

    #[error("Admixture node '{node}' does not have exactly 2 incoming edges")]
    InvalidAdmixtureNode { node: String },

    #[error("Cycle detected in graph at node '{node}'")]
    CycleDetected { node: String },

    // === Configuration Errors ===
    #[error("Invalid configuration: {details}")]
    InvalidConfig { details: String },

    #[error("Invalid thread count: {0} (must be 1 to 1024)")]
    InvalidThreadCount(usize),

    // === General Validation ===
    #[error("Validation error: {0}")]
    Validation(String),

    // === General Numerical ===
    #[error("Numerical error: {0}")]
    Numerical(String),

    // === Wrapped Errors ===
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Linear algebra error: {context}: {details}")]
    LinAlg { context: String, details: String },
}

impl PopmixError {
    /// Create an I/O error with file context.
    pub fn io_with_path(err: std::io::Error, path: &str) -> Self {
        match err.kind() {
            std::io::ErrorKind::NotFound => PopmixError::FileNotFound {
                path: path.to_string(),
            },
            _ => PopmixError::Io(err),
        }
    }

    /// Create a parse error with location context.
    pub fn parse(file: impl Into<String>, line: usize, details: impl Into<String>) -> Self {
        PopmixError::ParseError {
            file: file.into(),
            line,
            details: details.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = PopmixError::PopulationNotFound("French".to_string());
        assert_eq!(err.to_string(), "Population not found: 'French'");
    }

    #[test]
    fn test_error_from_io() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file not found");
        let err: PopmixError = io_err.into();
        assert!(matches!(err, PopmixError::Io(_)));
    }

    #[test]
    fn test_io_with_path() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "not found");
        let err = PopmixError::io_with_path(io_err, "/path/to/file.bed");
        assert!(matches!(err, PopmixError::FileNotFound { .. }));
    }
}
