//! Thread pool configuration for parallel computation.
//!
//! Popmix defaults to 4 threads for parallel operations.

use crate::error::{PopmixError, Result};

/// Default number of threads for parallel operations.
pub const DEFAULT_THREADS: usize = 4;

/// Initialize the global thread pool with the specified number of threads.
///
/// This should be called once at program startup. If not called, rayon will
/// use its default (number of logical CPUs).
///
/// # Arguments
///
/// * `threads` - Number of threads to use. If `None`, uses `DEFAULT_THREADS` (4).
///
/// # Returns
///
/// Returns `true` if the pool was initialized, `false` if it was already initialized.
/// Returns an error only if the thread count is 0.
///
/// # Example
///
/// ```
/// use popmix_core::init_thread_pool;
///
/// // Use default 4 threads
/// let initialized = init_thread_pool(None).unwrap();
///
/// // Or specify a custom count
/// // let initialized = init_thread_pool(Some(8)).unwrap();
/// ```
pub fn init_thread_pool(threads: Option<usize>) -> Result<bool> {
    let n = threads.unwrap_or(DEFAULT_THREADS);

    if n == 0 {
        return Err(PopmixError::InvalidThreadCount(0));
    }

    match rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()
    {
        Ok(()) => Ok(true),
        Err(_) => {
            // Pool already initialized - this is not an error
            Ok(false)
        }
    }
}

/// Get the current number of threads in the global thread pool.
///
/// This returns the actual number of threads rayon is using, which may differ
/// from the requested count if the pool was already initialized.
pub fn current_threads() -> usize {
    rayon::current_num_threads()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_threads() {
        assert_eq!(DEFAULT_THREADS, 4);
    }

    #[test]
    fn test_invalid_thread_count() {
        let result = init_thread_pool(Some(0));
        assert!(matches!(result, Err(PopmixError::InvalidThreadCount(0))));
    }

    #[test]
    fn test_current_threads() {
        // Just verify it returns a reasonable value
        let n = current_threads();
        assert!(n >= 1);
    }
}
