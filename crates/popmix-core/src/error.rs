//! Error types for popmix.
//!
//! All errors are explicit and propagated - NO exception swallowing. An error
//! about a file names the file; one about the caller's input says which input
//! and what was wrong with it, in a variant a caller can match on.

use std::path::{Path, PathBuf};

use thiserror::Error;

/// Result type alias for popmix operations.
pub type Result<T> = std::result::Result<T, PopmixError>;

/// Every error popmix reports.
///
/// New variants may be added in minor releases.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum PopmixError {
    // === Files ===
    #[error("File not found: {}", path.display())]
    FileNotFound { path: PathBuf },

    /// Reading or writing `path` failed; `source` is the operating system's
    /// error.
    #[error("{}: {source}", path.display())]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("Invalid file format in '{}': {details}", path.display())]
    InvalidFileFormat { path: PathBuf, details: String },

    #[error(
        "Unexpected file size for '{}': expected {expected} bytes, got {actual}",
        path.display()
    )]
    UnexpectedFileSize {
        path: PathBuf,
        expected: u64,
        actual: u64,
    },

    #[error("Failed to parse '{}' at line {line}: {details}", path.display())]
    ParseError {
        path: PathBuf,
        line: usize,
        details: String,
    },

    // === The caller's input ===
    #[error("Population not found: '{0}'")]
    PopulationNotFound(String),

    #[error("Duplicate population: '{0}'")]
    DuplicatePopulation(String),

    /// Fewer of something (populations, SNPs, blocks, samples) than the
    /// analysis needs.
    #[error("{context} needs at least {required} {what}, got {actual}")]
    TooFew {
        context: String,
        what: String,
        required: usize,
        actual: usize,
    },

    #[error("No SNPs remaining after filtering")]
    NoSnpsRemaining,

    #[error("{what} {index} out of range for {len}")]
    IndexOutOfRange {
        what: String,
        index: usize,
        len: usize,
    },

    #[error("Dimension mismatch in {context}: expected {expected}, got {actual}")]
    DimensionMismatch {
        context: String,
        expected: String,
        actual: String,
    },

    /// An argument outside what it may be (a bin size, a rank, a range).
    #[error("Invalid {name}: {details}")]
    InvalidParameter { name: String, details: String },

    #[error("Invalid thread count: {0} (must be 1 to 1024)")]
    InvalidThreadCount(usize),

    #[error("Invalid graph structure: {details}")]
    InvalidGraph { details: String },

    /// The data cannot support the analysis: no SNPs shared by two files,
    /// no block with a usable statistic, unsorted SNPs.
    #[error("{details}")]
    InvalidData { details: String },

    // === The computation ===
    #[error("Matrix is singular or near-singular in {context}")]
    SingularMatrix { context: String },

    #[error("Optimization failed to converge after {iterations} iterations in {context}")]
    ConvergenceFailure { iterations: usize, context: String },

    #[error("No feasible solution in {context}")]
    Infeasible { context: String },

    /// A numerical step failed: an SVD that did not converge, a non-finite
    /// value where a finite one is needed.
    #[error("{context}: {details}")]
    NumericalFailure { context: String, details: String },

    /// The thread pool could not be started.
    #[error("Thread pool: {0}")]
    ThreadPool(String),

    /// A bug in popmix, not a problem with the input.
    #[error("Internal error (please report it): {0}")]
    Internal(String),
}

/// What kind of failure an error is, for callers that map errors onto their
/// own (the Python bindings onto Python's exception types). A closed set:
/// a new [`PopmixError`] variant takes one of these kinds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorKind {
    /// A file that is not there.
    NotFound,
    /// Reading or writing a file failed.
    Io,
    /// The caller's input: a bad argument, file or population name.
    Input,
    /// The computation failed on valid input (a singular system).
    Computation,
    /// A bug in popmix.
    Internal,
}

impl PopmixError {
    /// What kind of failure this is.
    pub fn kind(&self) -> ErrorKind {
        use PopmixError as E;
        match self {
            E::FileNotFound { .. } => ErrorKind::NotFound,
            E::Io { source, .. } if source.kind() == std::io::ErrorKind::NotFound => {
                ErrorKind::NotFound
            }
            E::Io { .. } | E::ThreadPool(_) => ErrorKind::Io,
            E::InvalidFileFormat { .. }
            | E::UnexpectedFileSize { .. }
            | E::ParseError { .. }
            | E::PopulationNotFound(_)
            | E::DuplicatePopulation(_)
            | E::TooFew { .. }
            | E::NoSnpsRemaining
            | E::IndexOutOfRange { .. }
            | E::DimensionMismatch { .. }
            | E::InvalidParameter { .. }
            | E::InvalidThreadCount(_)
            | E::InvalidGraph { .. }
            | E::InvalidData { .. } => ErrorKind::Input,
            E::SingularMatrix { .. }
            | E::ConvergenceFailure { .. }
            | E::Infeasible { .. }
            | E::NumericalFailure { .. } => ErrorKind::Computation,
            E::Internal(_) => ErrorKind::Internal,
        }
    }

    /// An I/O error while reading or writing `path`: `FileNotFound` when the
    /// file is missing, otherwise `Io`, keeping the operating system's error.
    pub fn io(source: std::io::Error, path: impl AsRef<Path>) -> Self {
        let path = path.as_ref().to_path_buf();
        match source.kind() {
            std::io::ErrorKind::NotFound => PopmixError::FileNotFound { path },
            _ => PopmixError::Io { path, source },
        }
    }

    /// A parse error at `line` of `path`.
    pub fn parse(path: impl AsRef<Path>, line: usize, details: impl Into<String>) -> Self {
        PopmixError::ParseError {
            path: path.as_ref().to_path_buf(),
            line,
            details: details.into(),
        }
    }

    /// A malformed file.
    pub fn format(path: impl AsRef<Path>, details: impl Into<String>) -> Self {
        PopmixError::InvalidFileFormat {
            path: path.as_ref().to_path_buf(),
            details: details.into(),
        }
    }

    /// Too few of `what` for `context`.
    pub fn too_few(
        context: impl Into<String>,
        what: impl Into<String>,
        required: usize,
        actual: usize,
    ) -> Self {
        PopmixError::TooFew {
            context: context.into(),
            what: what.into(),
            required,
            actual,
        }
    }

    /// An argument `name` that is not allowed, and why.
    pub fn parameter(name: impl Into<String>, details: impl Into<String>) -> Self {
        PopmixError::InvalidParameter {
            name: name.into(),
            details: details.into(),
        }
    }

    /// Data that cannot support the analysis.
    pub fn data(details: impl Into<String>) -> Self {
        PopmixError::InvalidData {
            details: details.into(),
        }
    }

    /// A failed numerical step in `context`.
    pub fn numerical(context: impl Into<String>, details: impl Into<String>) -> Self {
        PopmixError::NumericalFailure {
            context: context.into(),
            details: details.into(),
        }
    }

    /// `index` of `what` past the `len` there are.
    pub fn out_of_range(what: impl Into<String>, index: usize, len: usize) -> Self {
        PopmixError::IndexOutOfRange {
            what: what.into(),
            index,
            len,
        }
    }
}

/// Names the file in an I/O error: `File::open(p).with_path(p)?`.
pub trait IoResultExt<T> {
    /// The error as [`PopmixError::io`] for `path`.
    fn with_path(self, path: impl AsRef<Path>) -> Result<T>;
}

impl<T> IoResultExt<T> for std::result::Result<T, std::io::Error> {
    fn with_path(self, path: impl AsRef<Path>) -> Result<T> {
        self.map_err(|e| PopmixError::io(e, path))
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error as _;

    use super::*;

    #[test]
    fn a_population_error_names_the_population() {
        let err = PopmixError::PopulationNotFound("French".to_string());
        assert_eq!(err.to_string(), "Population not found: 'French'");
    }

    #[test]
    fn a_missing_file_is_file_not_found_with_its_path() {
        let io = std::io::Error::new(std::io::ErrorKind::NotFound, "no such file");
        let err = PopmixError::io(io, "/data/x.bed");
        assert!(
            matches!(&err, PopmixError::FileNotFound { path } if path == Path::new("/data/x.bed"))
        );
        assert_eq!(err.to_string(), "File not found: /data/x.bed");
    }

    #[test]
    fn other_io_errors_keep_the_path_and_the_cause() {
        let io = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "permission denied");
        let err = PopmixError::io(io, "/data/x.bed");
        assert_eq!(err.to_string(), "/data/x.bed: permission denied");
        let source = err.source().expect("the io error is the source");
        assert_eq!(source.to_string(), "permission denied");
        assert!(
            matches!(err, PopmixError::Io { ref path, .. } if path == Path::new("/data/x.bed"))
        );
    }

    #[test]
    fn with_path_names_the_file() {
        let read = std::fs::read("/nonexistent/popmix-test").with_path("/nonexistent/popmix-test");
        assert_eq!(
            read.unwrap_err().to_string(),
            "File not found: /nonexistent/popmix-test"
        );
    }

    #[test]
    fn each_error_has_a_kind() {
        let gone = std::io::Error::new(std::io::ErrorKind::NotFound, "gone");
        let denied = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "denied");
        assert_eq!(PopmixError::io(gone, "a").kind(), ErrorKind::NotFound);
        assert_eq!(PopmixError::io(denied, "a").kind(), ErrorKind::Io);
        assert_eq!(
            PopmixError::PopulationNotFound("X".into()).kind(),
            ErrorKind::Input
        );
        assert_eq!(
            PopmixError::too_few("qpAdm", "right populations", 2, 1).kind(),
            ErrorKind::Input
        );
        assert_eq!(
            PopmixError::numerical("qpAdm", "SVD failed").kind(),
            ErrorKind::Computation
        );
        assert_eq!(
            PopmixError::Internal("x".into()).kind(),
            ErrorKind::Internal
        );
    }

    #[test]
    fn input_errors_say_what_was_wrong() {
        assert_eq!(
            PopmixError::too_few("qpAdm", "right populations", 2, 1).to_string(),
            "qpAdm needs at least 2 right populations, got 1"
        );
        assert_eq!(
            PopmixError::out_of_range("sample", 7, 5).to_string(),
            "sample 7 out of range for 5"
        );
        assert_eq!(
            PopmixError::parameter("beam_width", "must be at least 1").to_string(),
            "Invalid beam_width: must be at least 1"
        );
        assert_eq!(
            PopmixError::parse("a.snp", 3, "bad position").to_string(),
            "Failed to parse 'a.snp' at line 3: bad position"
        );
    }
}
